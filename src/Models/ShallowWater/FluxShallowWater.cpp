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

#include "FluxShallowWater.h"

//***********************************************************************

FluxShallowWater::FluxShallowWater() : m_mass(0.), m_momentum(0.) {}

//***********************************************************************

FluxShallowWater::~FluxShallowWater() {}

//***********************************************************************

void FluxShallowWater::addFlux(double coefA)
{
  m_mass     += coefA * static_cast<FluxShallowWater*>(fluxBuff)->m_mass;
  m_momentum += coefA * static_cast<FluxShallowWater*>(fluxBuff)->m_momentum;
}

//***********************************************************************

void FluxShallowWater::addFlux(Flux* flux)
{
  m_mass     += static_cast<FluxShallowWater*>(flux)->m_mass;
  m_momentum += static_cast<FluxShallowWater*>(flux)->m_momentum;
}

//***********************************************************************

void FluxShallowWater::subtractFlux(double coefA)
{
  m_mass     -= coefA * static_cast<FluxShallowWater*>(fluxBuff)->m_mass;
  m_momentum -= coefA * static_cast<FluxShallowWater*>(fluxBuff)->m_momentum;
}

//***********************************************************************

void FluxShallowWater::multiply(double scalar)
{
  m_mass     *= scalar;
  m_momentum *= scalar;
}

//***********************************************************************

void FluxShallowWater::setBufferFlux(Cell& cell) { static_cast<FluxShallowWater*>(fluxBuff)->buildCons(cell.getPhases(), cell.getMixture()); }

//***********************************************************************

void FluxShallowWater::buildCons(Phase** phases, Mixture* /*mixture*/)
{
  Phase* phase{phases[0]};

  m_mass     = phase->getHeight();
  m_momentum = m_mass * phase->getVelocity();
}

//***********************************************************************

void FluxShallowWater::buildPrim(Phase** phases, Mixture* /*mixture*/)
{
  Phase* phase{phases[0]};
  double pressure{0.}, soundSpeed{0.};
  Eos* eos{phase->getEos()};

  phase->setHeight(m_mass);
  phase->setVelocity(m_momentum.getX() / m_mass, m_momentum.getY() / m_mass, 0.);
  //Erasing small velocity variations
  if (std::fabs(phase->getU()) < 1.e-8) phase->setU(0.);
  if (std::fabs(phase->getV()) < 1.e-8) phase->setV(0.);

  pressure = eos->computePressure(m_mass);
  phase->setPressure(pressure);
  soundSpeed = eos->computeSoundSpeed(m_mass);
  phase->setSoundSpeed(soundSpeed);
}

//***********************************************************************

void FluxShallowWater::setToZero()
{
  m_mass     = 0.;
  m_momentum = 0.;
}

//***********************************************************************
