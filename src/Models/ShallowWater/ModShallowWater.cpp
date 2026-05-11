//
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-.
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| |
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | |
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  |
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)|
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_)
//      (__)              (_)      (__)     (__)     (__)
//      Official webSite: https://code-mphi.github.io/ECOGEN/
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names
//  are listed in the copyright file included with this source
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published
//  by the Free Software Foundation, either version 3 of the License,
//  or (at your option) any later version.
//
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).
//  If not, see <http://www.gnu.org/licenses/>.

#include "ModShallowWater.h"

const std::string ModShallowWater::NAME = "SHALLOWWATER";

//****************************************************************************

ModShallowWater::ModShallowWater(const int& numbTransports) : Model(NAME, numbTransports)
{
  fluxBuff = new FluxShallowWater();
  for (int i = 0; i < 3; i++) {
    sourceCons.push_back(new FluxShallowWater());
  }
}

//****************************************************************************

ModShallowWater::~ModShallowWater()
{
  delete fluxBuff;
  for (int i = 0; i < 3; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//****************************************************************************

void ModShallowWater::allocateCons(Flux** cons) { *cons = new FluxShallowWater; }

//***********************************************************************

void ModShallowWater::allocatePhase(Phase** phase) { *phase = new PhaseShallowWater; }

//***********************************************************************

void ModShallowWater::allocateMixture(Mixture** mixture) { *mixture = new MixShallowWater; }

//***********************************************************************

void ModShallowWater::fulfillState(Phase** phases, Mixture* /*mixture*/)
{
  // Compute pressure and sound speed from EoS (velocity is not used but required for sake of genericity).
  phases[0]->extendedCalculusPhase(phases[0]->getVelocity());
}

//****************************************************************************
//********************* Cell to cell Riemann solvers *************************
//****************************************************************************

void ModShallowWater::solveRiemannIntern(
  Cell& cellLeft, Cell& cellRight, const double& dxLeft, const double& dxRight, double& dtMax, std::vector<double>& /*boundData*/) const
{
  double cL, cR, sL, sR;
  double uL, uR, vL, vR, hL, hR, pL, pR;

  Phase* phaseLeft{cellLeft.getPhase(0)};   // Get phase '0' (our only phase) of left cell
  Phase* phaseRight{cellRight.getPhase(0)}; // Get phase '0' (our only phase) of right cell

  // Get left state
  uL = phaseLeft->getU();
  vL = phaseLeft->getV();
  hL = phaseLeft->getHeight();
  pL = phaseLeft->getPressure();
  cL = phaseLeft->getSoundSpeed();

  // Get right state
  uR = phaseRight->getU();
  vR = phaseRight->getV();
  hR = phaseRight->getHeight();
  pR = phaseRight->getPressure();
  cR = phaseRight->getSoundSpeed();

  sL = std::min(uL - cL, uR - cR);
  sR = std::max(uR + cR, uL + cL);

  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));
  if (std::fabs(sR) > 1.e-3) dtMax = std::min(dtMax, dxRight / std::fabs(sR));

  //compute left and right mass flow rates and sM
  double mL(hL * (sL - uL)), mR(hR * (sR - uR));
  double sM((pR - pL + mL * uL - mR * uR) / (mL - mR));
  if (std::fabs(sM) < 1.e-8) sM = 0.;

  if (sL > 0.) {
    static_cast<FluxShallowWater*>(fluxBuff)->m_mass = hL * uL;
    static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setX(hL * uL * uL + pL);
    static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY(hL * vL * uL);
  }
  else if (sR < 0.) {
    static_cast<FluxShallowWater*>(fluxBuff)->m_mass = hR * uR;
    static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setX(hR * uR * uR + pR);
    static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY(hR * vR * uR);
  }
  // HLLC
  else if (std::fabs(sR - sL) > 1.e-3) {
#warning to understand
    static_cast<FluxShallowWater*>(fluxBuff)->m_mass = (hR * uR * sL - hL * uL * sR + sL * sR * (hL - hR)) / (sL - sR);
    double momFluxX = ((hR * uR * uR + pR) * sL - (hL * uL * uL + pL) * sR + sL * sR * (hL * uL - hR * uR)) / (sL - sR);
    //Correction for W-P scheme
    //momFluxX += 0.5*phaseLeft->getEos()->getG()*(hL+hR)*(zR-zL);
    static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setX(momFluxX);
    if (sM >= 0) {
      static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY(static_cast<FluxShallowWater*>(fluxBuff)->m_mass * vL);
    }
    else {
      static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY(static_cast<FluxShallowWater*>(fluxBuff)->m_mass * vR);
    }
  }
  else {
    std::cout << "Unexpected, sL = sR = " << sL << std::endl;
    exit(EXIT_FAILURE);
  }

  // HLLC Euler
  // else if (sM >= 0.) {
  //   double pStar   = mL * (sM - uL) + pL;
  //   double rhoStar = mL / (sL - sM);

  //   static_cast<FluxShallowWater*>(fluxBuff)->m_mass = rhoStar * sM;
  //   static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setX(rhoStar * sM * sM + pStar);
  //   static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY(rhoStar * sM * vL);

  // }
  // else {
  //   double pStar   = mR * (sM - uR) + pR;
  //   double rhoStar = mR / (sR - sM);
  //   double Estar   = ER + (sM - uR) * (sM + pR / mR);

  //   static_cast<FluxShallowWater*>(fluxBuff)->m_mass = rhoStar * sM;
  //   static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setX(rhoStar * sM * sM + pStar);
  //   static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY(rhoStar * sM * vR);
  //   static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setZ(rhoStar * sM * wR);
  //   static_cast<FluxShallowWater*>(fluxBuff)->m_energ = (rhoStar * Estar + pStar) * sM;

  //   // Boundary data for output
  //   boundData[VarBoundary::p]    = pStar;
  //   boundData[VarBoundary::rho]  = rhoStar;
  //   boundData[VarBoundary::velU] = sM;
  //   boundData[VarBoundary::velV] = vR;
  //   boundData[VarBoundary::velW] = wR;
  // }

  //Contact discontinuity
  static_cast<FluxShallowWater*>(fluxBuff)->m_sM = sM;
  static_cast<FluxShallowWater*>(fluxBuff)->m_hL = hL;
  static_cast<FluxShallowWater*>(fluxBuff)->m_hR = hR;
}

//****************************************************************************
//************** Half Riemann solvers for boundary conditions** **************
//****************************************************************************

void ModShallowWater::solveRiemannWall(Cell& cellLeft, const double& dxLeft, double& dtMax, std::vector<double>& /*boundData*/) const
{
  double cL, sL;
  double uL, pL, hL;
  double pStar{0.};

  // Get left state
  Phase* phaseLeft{cellLeft.getPhase(0)};
  uL = phaseLeft->getU();
  pL = phaseLeft->getPressure();
  hL = phaseLeft->getHeight();
  cL = phaseLeft->getSoundSpeed();

  sL = std::min(uL - cL, -uL - cL);
  if (std::fabs(sL) > 1.e-3) dtMax = std::min(dtMax, dxLeft / std::fabs(sL));

  pStar = hL * uL * (uL - sL) + pL;

  static_cast<FluxShallowWater*>(fluxBuff)->m_mass = 0.;
  static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setX(pStar);
  static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.setY(0.);

  //Contact discontinuity velocity
  static_cast<FluxShallowWater*>(fluxBuff)->m_sM = 0.;
  static_cast<FluxShallowWater*>(fluxBuff)->m_hL = hL;
}

//****************************************************************************
//******************************* Accessors **********************************
//****************************************************************************

double ModShallowWater::selectScalar(Phase** phases, Mixture* /*mixture*/, Transport* /*transports*/, Variable nameVariable, int /*num*/) const
{
  switch (nameVariable) {
  case Variable::pressure:
    return phases[0]->getPressure();
    break;
  case Variable::height:
    return phases[0]->getHeight();
    break;
  case Variable::velocityU:
    return phases[0]->getU();
    break;
  case Variable::velocityV:
    return phases[0]->getV();
    break;
  case Variable::velocityMag:
    return phases[0]->getVelocity().norm();
    break;
  default:
    Errors::errorMessage("nameVariable unknown in selectScalar.");
    return 0;
    break;
  }
}

//****************************************************************************
//***************************** Others methods *******************************
//****************************************************************************

void ModShallowWater::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  static_cast<FluxShallowWater*>(fluxBuff)->m_momentum.reverseProjection(normal, tangent, binormal);
}

//****************************************************************************
