/* -------------------------------------------------------------------------- *
 *                    OpenSim:  ThreeCompartmentFatigue.cpp                   *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "ThreeCompartmentFatigue.h"

using namespace OpenSim;

static const char* ACTIVE_NAME = "active_fraction";
static const char* FATIGUED_NAME = "fatigued_fraction";

ThreeCompartmentFatigue::ThreeCompartmentFatigue()
{
    constructProperties();
}

void ThreeCompartmentFatigue::constructProperties()
{
    constructProperty_fatigue_rate(0.00755);
    constructProperty_recovery_rate(0.00075);
    constructProperty_rest_recovery_multiplier(15.0);
    constructProperty_development_factor(10.0);
    constructProperty_relaxation_factor(10.0);
    constructProperty_rest_target_load(0.0);
    constructProperty_initial_active_fraction(0.0);
    constructProperty_initial_fatigued_fraction(0.0);
}

double ThreeCompartmentFatigue::getActiveFraction(const SimTK::State& s) const
{
    return getStateVariableValue(s, ACTIVE_NAME);
}

double ThreeCompartmentFatigue::getFatiguedFraction(const SimTK::State& s) const
{
    return getStateVariableValue(s, FATIGUED_NAME);
}

double ThreeCompartmentFatigue::getRestingFraction(const SimTK::State& s) const
{
    return 1.0 - getActiveFraction(s) - getFatiguedFraction(s);
}

double ThreeCompartmentFatigue::getForceCapacity(const SimTK::State& s) const
{
    return 1.0 - getFatiguedFraction(s);
}

void ThreeCompartmentFatigue::setActiveFraction(
        SimTK::State& s, double value) const
{
    setStateVariableValue(s, ACTIVE_NAME, value);
}

void ThreeCompartmentFatigue::setFatiguedFraction(
        SimTK::State& s, double value) const
{
    setStateVariableValue(s, FATIGUED_NAME, value);
}

void ThreeCompartmentFatigue::calcDerivatives(double ma, double mf,
        double targetLoad, double& activeRate, double& fatiguedRate) const
{
    const double mr = 1.0 - ma - mf;
    const double need = targetLoad - ma;

    double c;
    if (ma < targetLoad) {
        c = (mr > need) ? get_development_factor() * need
                        : get_development_factor() * mr;
    } else {
        c = get_relaxation_factor() * need;
    }

    const double recovery = (targetLoad <= get_rest_target_load())
            ? get_recovery_rate() * get_rest_recovery_multiplier()
            : get_recovery_rate();

    activeRate = c - get_fatigue_rate() * ma;
    fatiguedRate = get_fatigue_rate() * ma - recovery * mf;
}

void ThreeCompartmentFatigue::extendAddToSystem(
        SimTK::MultibodySystem& system) const
{
    Super::extendAddToSystem(system);
    addStateVariable(ACTIVE_NAME);
    addStateVariable(FATIGUED_NAME);
}

void ThreeCompartmentFatigue::extendInitStateFromProperties(
        SimTK::State& s) const
{
    Super::extendInitStateFromProperties(s);
    setActiveFraction(s, get_initial_active_fraction());
    setFatiguedFraction(s, get_initial_fatigued_fraction());
}

void ThreeCompartmentFatigue::extendSetPropertiesFromState(
        const SimTK::State& s)
{
    Super::extendSetPropertiesFromState(s);
    set_initial_active_fraction(getActiveFraction(s));
    set_initial_fatigued_fraction(getFatiguedFraction(s));
}

void ThreeCompartmentFatigue::computeStateVariableDerivatives(
        const SimTK::State& s) const
{
    const Muscle& muscle = getConnectee<Muscle>("muscle");
    const double ma = getActiveFraction(s);
    const double mf = getFatiguedFraction(s);
    // Target load is the activation as a fraction of the unfatigued force.
    const double targetLoad = muscle.getActivation(s) * (1.0 - mf);

    double activeRate, fatiguedRate;
    calcDerivatives(ma, mf, targetLoad, activeRate, fatiguedRate);

    setStateVariableDerivativeValue(s, ACTIVE_NAME, activeRate);
    setStateVariableDerivativeValue(s, FATIGUED_NAME, fatiguedRate);
}
