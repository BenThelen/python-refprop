#-------------------------------------------------------------------------------
#Name:            multiRP
#Purpose:         allowing multiple refprop calculations simultaneously with
#                 multiple CPU support using python multiprocessing module
#
#Author:          Thelen, B.J.
#                 thelen_ben@yahoo.com
#-------------------------------------------------------------------------------

'''Multiprocessing can be done with refprop.py module provided that the setup
details of each refprop routine called simultaneously is identical.

This module manage and control call of various routines of different setup
details such that only routines of identical setup details are being called
simultanously. The no. cores and threads of cpu matches the the maximum multiple
refprop calculations to be performed simultanously.

The aim of this module is to gain time of fluid cycle calculations.

All the same functions as found in the refprop module are available with
additional input of the setup details and porting of the multiprocessing
variables

Note for windows users, initialtion of multirefprop does require some time
which could render the time gain. On my system the initiation difference ratio
between Linux and Windows is 14.04365604 times (in favour for Linux)'''
import os
import sys
import refprop
import time
import multiprocessing as mp
from decimal import Decimal

#input declarations
RefpropError = refprop.RefpropError
RefpropdllError = refprop.RefpropdllError
RefpropicompError = refprop.RefpropicompError
RefpropinputError = refprop.RefpropinputError
RefpropnormalizeError = refprop.RefpropnormalizeError
RefproproutineError = refprop.RefproproutineError
RefpropWarning = refprop.RefpropWarning
RefpropdllWarning = refprop.RefpropdllWarning

#Classes
class _MultiRefProp(mp.Process):
    '''initiate multiprocessing for refprop,

    this function needs to be called prior to usage of multiRP and under
    "if __name__ == '__main__':" enable to function properly'''
    def __init__(self):
        self.process = mp.Process
        self.mgr = mp.Manager()
        self.sem = mp.Semaphore(mp.cpu_count())
        self.result = self.mgr.dict() #result dict
        self.ppipe, self.cpipe = mp.Pipe() #parent pipe, #child pipe

class MultiRPError(RefpropError):
    'General error for multiRP'
    pass

class SetupError(MultiRPError):
    'Raise input error when Setups are blocked'
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class MultiRPInputError(MultiRPError):
    'Raise input error when Setups are blocked'
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class MultiRPChildError(MultiRPError):
    'Raise input error when Setups are blocked'
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


#Classes from refprop.py
class FluidModel(refprop.FluidModel):
    '''return array of fluid model'''
    pass

class SetWarning():
    'Return RefpropdllWarning status (on / off)'
    def __repr__(self):
        return str(refprop.SetWarning())
    @staticmethod
    def on(prop=None, mRP=None):
        'Sets RefpropdllWarning on, initiate Error on Refpropdll ierr value < 0'
        global _inhibt_setuperror
        _inhibt_setuperror = True
        def _rpfunc():
            global _inhibt_setuperror, _setupprop
            _inhibt_setuperror = False
            prop = refprop.SetWarning.on()
            #adjust _setupprop
            _setupprop = setup_details(prop)
            return prop
        return _rpfunc_handler(prop, mRP, _rpfunc)
    @staticmethod
    def off(prop=None, mRP=None):
        'Sets RefpropdllWarning off, no Error raised on Refpropdll ierr value < 0'
        global _inhibt_setuperror
        _inhibt_setuperror = True
        def _rpfunc():
            global _inhibt_setuperror, _setupprop
            _inhibt_setuperror = False
            prop = refprop.SetWarning.off()
            #adjust _setupprop
            _setupprop = setup_details(prop)
            return prop
        return _rpfunc_handler(prop, mRP, _rpfunc)

class SetError():
    'Return RefpropdllError status (on / off)'
    def __repr__(self):
        return str(refprop.SetError())
    @staticmethod
    def on(prop=None, mRP=None):
        'Sets RefpropdllError on, initiate Error on Refpropdll ierr value != 0'
        global _inhibt_setuperror
        _inhibt_setuperror = True
        def _rpfunc():
            global _inhibt_setuperror, _setupprop
            _inhibt_setuperror = False
            prop = refprop.SetError.on()
            #adjust _setupprop
            _setupprop = setup_details(prop)
            return prop
        return _rpfunc_handler(prop, mRP, _rpfunc)
    @staticmethod
    def off(prop=None, mRP=None):
        'Sets RefpropdllError off, no Error raised on Refpropdll ierr value != 0'
        global _inhibt_setuperror
        _inhibt_setuperror = True
        def _rpfunc():
            global _inhibt_setuperror, _setupprop
            _inhibt_setuperror = False
            prop = refprop.SetError.off()
            #adjust _setupprop
            _setupprop = setup_details(prop)
            return prop
        return _rpfunc_handler(prop, mRP, _rpfunc)

class SetErrorDebug():
    'Return SetErrorDebug status (on / off)'
    def __repr__(self):
        return str(refprop.SetErrorDebug())
    @staticmethod
    def on(prop=None, mRP=None):
        'Sets error debug mode on, displays error message only'
        global _inhibt_setuperror
        _inhibt_setuperror = True
        def _rpfunc():
            global _inhibt_setuperror, _setupprop
            _inhibt_setuperror = False
            prop = refprop.SetErrorDebug.on()
            #adjust _setupprop
            _setupprop = setup_details(prop)
            return prop
        return _rpfunc_handler(prop, mRP, _rpfunc)
    @staticmethod
    def off(prop=None, mRP=None):
        'Sets error debug mode off, displays error message only'
        global _inhibt_setuperror
        _inhibt_setuperror = True
        def _rpfunc():
            global _inhibt_setuperror, _setupprop
            _inhibt_setuperror = False
            prop = refprop.SetErrorDebug.off()
            #adjust _setupprop
            _setupprop = setup_details(prop)
            return prop
        return _rpfunc_handler(prop, mRP, _rpfunc)


#functions
def multirefprop():
    'initiate multiprocessing variables'
    global mRP
    #identify parent process (only parent can call this)
    if mp.current_process()._parent_pid == None:
        #only call if never called
        if not 'mRP' in globals():
            _mRP = _MultiRefProp()
            mRP = {'sem':_mRP.sem,
                   'process':_mRP.process,
                   'result':_mRP.result, 'ppipe':_mRP.ppipe,
                   'cpipe':_mRP.cpipe}
    else:
        raise MultiRPChildError(
        'Only parent process can initiate mRP variables')

    return mRP

def run_mRP(processlist):
    'start, close and check child processes'
    #start child processes
    for each in processlist:
        each.start()
    #close child processes
    for each in processlist:
        each.join()
    #check for errors in child processes
    for each in processlist:
        if each.exitcode != 0:
            raise MultiRPChildError('child error in ' + each.name +
                                    ' exitcode = ' + str(each.exitcode))


def _multirefprop(_mRP):
    'Set parent multiprocessing variables as globals'
    global sem, process, result, mRP, ppipe, cpipe
    sem = _mRP['sem']
    result = _mRP['result']
    process = _mRP['process']
    mRP = _mRP
    cpipe = _mRP['cpipe']
    ppipe = _mRP['ppipe']

def ppip():
    'return parent pipe'
    if 'mRP' in globals():
        return mRP['pipe']
    else:
        raise MultiRPInputError('parent pipe can not be return without a' +
                                 ' multirefprop() call')

def _rpfunc_handler(prop, mRP, _rpfunc):
    global _setupprop
    #verify either prop has value or value has been set before'
    if prop == None \
    and ('_setupprop' not in dir(refprop) \
    or refprop._setupprop == {}) \
    and '_inhibt_setuperror' in globals() \
    and not _inhibt_setuperror:
        raise MultiRPInputError(
            'First refprop function in child needs setup' +
            ' details in "prop" to be defined' +
            ' or parent process needs setup input first')
    if mRP != None:
        #set mRP values as globals
        _multirefprop(mRP)
    if prop != None:
        #store prop as global for continue function call in same child
        _setupprop = setup_details(prop)
    #identify child process
    if mp.current_process()._parent_pid != None:
        if prop == None:
            if '_setupprop' not in globals() or _setupprop == None:
                raise MultiRPInputError(
                'child process needs "prop" input for first refprop function')
            #function call to inherite prop values from previous call
            prop = _setupprop
        #raise error if child process is initiated without mRP settings
        if not 'mRP' in globals():
            raise MultiRPInputError(
                'child process needs "mRP" input for first refprop function')
        #resetup fluid
        refprop.resetup(prop)
        with sem:
            #run function, this is where the refprop function is called
            prps = _rpfunc()
        #return result in result.dict (for child processes only)
        result[process().name.rsplit(':', 1)[0]] = prps
        return prps
    #identify parent process
    elif mp.current_process()._parent_pid == None:
        #confirm active children
        if len(mp.active_children()) > 1:
            #raise error if children are active
            raise MultiRPInputError(
                'parent process can not proceed with ' +
                'child processes being handled. Wait untill ' +
                'child processes have joined. Use run_mRP ' +
                'command to start and join child processes.')
        #resetup
        if prop != None:
            refprop.resetup(prop)
        #run function
        return _rpfunc()

def _checksetupblock(name): ##################################################################perhaps this whole function can be ommitted, but carefully as settings are passed through
    if mp.current_process()._parent_pid != None:
        #raise error if child process request for setup
        raise SetupError('function "' + str(name) +
                          '" is blocked for multiRP child processes')
    elif len(mp.active_children()) > 1:
        #raise error if parent process request for setup while child process(es)
        #are active.
        raise SetupError(
            'function "' + str(name) +
            '" is blocked while multiRP child processes are active')


#functions from refprop.py
def setpath(path='c:/program files/refprop/'):
    '''Set Directory to refprop root containing fluids and mixtures.'''
    refprop.setpath(path)

def fluidlib():
    '''Displays all fluids and mixtures available on root directory.'''
    refprop.fluidlib()

def normalize(x):
    '''Normalize the sum of list x value's to 1'''
    return refprop.normalize(x)

def resetup(prop, force=False, mRP=None):
    """Resetup models and re-initialize arrays.
    Initiate child process with prop and mRP transfer"""
    #resetup only to be run in parent process, in child process the passed prop
    #value will be used for resetup
    def _rpfunc():
        return refprop.resetup(prop, force)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def setup_details(prop):
    '''Returns basic setup details.'''
    return refprop.setup_details(prop)

def getphase(prop, mRP=None):
    '''Return fluid phase'''
    def _rpfunc():
        return refprop.getphase(prop)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def setup_setting(mRP=None):
    """Returns current setup settings"""
    #identify parent process and confirm no active children
    if (mp.current_process()._parent_pid == None and
    len(mp.active_children()) <= 1):
        return refprop.setup_setting()
    #else multiprocessing active
    else:
        #check if mRP is included
        if 'mRP' not in globals():
            raise MultiRPInputError(
            '"mRP" input required when multiprocessing')
        #check if setup_prop is already populated (done through first call)
        if '_setupprop' not in globals() or _setupprop == None:
                raise MultiRPInputError(
                'child process needs "prop" input for first refprop function')
        #return _Setupprop (which contains setup-properties).
        return _setupprop

def _test():
    """execute detailed test run of multiRP"""
    import rptest
    rptest.settest('multiRP')

def test(criteria = 0.00001):
    '''verify that the user's computer is returning proper calculations'''
    global testresult
    refproptest = refprop.test(criteria)
    testresult = refprop.testresult
    return refproptest

def psliq(p, s, x, prop=None, mRP=None):
    """flsh1 calculations with boundery check, raise RefpropinputError when
    input is outside bounderies"""
    def _rpfunc():
        return refprop.psliq(p, s, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def psvap(p, s, x, prop=None, mRP=None):
    """flsh1 calculations with boundery check, raise RefpropinputError when
    input is outside bounderies"""
    def _rpfunc():
        return refprop.psvap(p, s, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def ps2ph(p, s, x, prop=None, mRP=None):
    """flsh2 calculations with boundery check, raise RefpropinputError when
    input is outside bounderies"""
    def _rpfunc():
        return refprop.ps2ph(p, s, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def phliq(p, h, x, prop=None, mRP=None):
    """flsh1 calculations with boundery check, raise RefpropinputError when
    input is outside bounderies"""
    def _rpfunc():
        return refprop.phliq(p, h, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def phvap(p, h, x, prop=None, mRP=None):
    """flsh1 calculations with boundery check, raise RefpropinputError when
    input is outside bounderies"""
    def _rpfunc():
        return refprop.phvap(p, h, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def ph2ph(p, h, x, prop=None, mRP=None):
    """flsh2 calculations with boundery check, raise RefpropinputError when
    input is outside bounderies"""
    def _rpfunc():
        return refprop.ph2ph(p, h, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def setup(hrf, *hfld, hfmix='HMX.BNC'):
    '''Define models and initialize arrays.'''
    _checksetupblock('setup')
    return refprop.setup(hrf, *hfld, hfmix=hfmix)

def setmod(htype, hmix, *hcomp):
    '''Set model(s) other than the NIST-recommended ('NBS') ones.'''
    _checksetupblock('setmod')
    refprop.setmod(htype, hmix, *hcomp)

def gerg04(ixflag=0):
    '''set the pure model(s) to those used by the GERG 2004 formulation.'''
    _checksetupblock('gerg04')
    refprop.gerg04(ixflag)

def setref(hrf='DEF', ixflag=1, x0=[1], h0=0, s0=0, t0=273, p0=100):
    '''set reference state enthalpy and entropy'''
    _checksetupblock('setref')
    return refprop.setref(hrf, ixflag, x0, h0, s0, t0, p0)

def critp(x, prop=None, mRP=None):
    '''critical parameters as a function of composition'''
    def _rpfunc():
        return refprop.critp(x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def therm(t, D, x, prop=None, mRP=None):
    '''Compute thermal quantities as a function of temperature, density and
    compositions using core functions'''
    def _rpfunc():
        return refprop.therm(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def therm0(t, D, x, prop=None, mRP=None):
    '''Compute ideal gas thermal quantities as a function of temperature,
    density and compositions using core functions.'''
    def _rpfunc():
        return refprop.therm0(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def residual (t, D, x, prop=None, mRP=None):
    '''compute the residual quantities as a function of temperature, density,
    and compositions (where the residual is the property minus the ideal gas
    portion)'''
    def _rpfunc():
        return refprop.residual (t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def therm2(t, D, x, prop=None, mRP=None):
    '''Compute thermal quantities as a function of temperature, density and
    compositions using core functions'''
    def _rpfunc():
        return refprop.therm2(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def therm3(t, D, x, prop=None, mRP=None):
    '''Compute miscellaneous thermodynamic properties'''
    def _rpfunc():
        return refprop.therm3(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def purefld(icomp=0):
    '''Change the standard mixture setup so that the properties of one fluid
    can be calculated as if SETUP had been called for a pure fluid.'''
    _checksetupblock('purefld')
    return refprop.purefld(icomp)

def name(icomp=1, prop=None, mRP=None):
    '''Provides name information for specified component'''
    def _rpfunc():
        return refprop.name(icomp)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def entro(t, D, x, prop=None, mRP=None):
    '''Compute entropy as a function of temperature, density and composition
    using core functions'''
    def _rpfunc():
        return refprop.entro(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def enthal(t, D, x, prop=None, mRP=None):
    '''Compute enthalpy as a function of temperature, density, and
    composition using core functions'''
    def _rpfunc():
        return refprop.enthal(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def cvcp(t, D, x, prop=None, mRP=None):
    '''Compute isochoric (constant volume) and isobaric (constant pressure)
    heat capacity as functions of temperature, density, and composition
    using core functions'''
    def _rpfunc():
       return refprop.cvcp(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def cvcpk(icomp, t, D, prop=None, mRP=None):
    '''Compute isochoric (constant volume) and isobaric (constant pressure)
    heat capacity as functions of temperature for a given component.'''
    def _rpfunc():
        return refprop.cvcpk(icomp, t, D)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def gibbs(t, D, x, prop=None, mRP=None):
    '''Compute residual Helmholtz and Gibbs free energy as a function of
    temperature, density, and composition using core functions'''
    def _rpfunc():
        return refprop.gibbs(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def ag(t, D, x, prop=None, mRP=None):
    '''Ccompute Helmholtz and Gibbs energies as a function of temperature,
    density, and composition.'''
    def _rpfunc():
        return refprop.ag(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def press(t, D, x, prop=None, mRP=None):
    '''Compute pressure as a function of temperature, density, and
    composition using core functions'''
    def _rpfunc():
        return refprop.press(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dpdd(t, D, x, prop=None, mRP=None):
    '''Compute partial derivative of pressure w.r.t. density at constant
    temperature as a function of temperature, density, and composition'''
    def _rpfunc():
        return refprop.dpdd(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dpddk(icomp, t, D, prop=None, mRP=None):
    '''Compute partial derivative of pressure w.r.t. density at constant
    temperature as a function of temperature and density for a specified
    component'''
    def _rpfunc():
        return refprop.dpddk(icomp, t, D)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dpdd2(t, D, x, prop=None, mRP=None):
    '''Compute second partial derivative of pressure w.r.t. density at
    const. temperature as a function of temperature, density, and
    composition.'''
    def _rpfunc():
        return refprop.dpdd2(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dpdt(t, D, x, prop=None, mRP=None):
    '''Compute partial derivative of pressure w.r.t. temperature at constant
    density as a function of temperature, density, and composition.'''
    def _rpfunc():
        return refprop.dpdt(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dpdtk(icomp, t, D, prop=None, mRP=None):
    '''Compute partial derivative of pressure w.r.t. temperature at constant
    density as a function of temperature and density for a specified
    component'''
    def _rpfunc():
        return refprop.dpdtk(icomp, t, D)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dddp(t, D, x, prop=None, mRP=None):
    '''ompute partial derivative of density w.r.t. pressure at constant
    temperature as a function of temperature, density, and composition.'''
    def _rpfunc():
        return refprop.dddp(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dddt(t, D, x, prop=None, mRP=None):
    '''Compute partial derivative of density w.r.t. temperature at constant
    pressure as a function of temperature, density, and composition.'''
    def _rpfunc():
        return refprop.dddt(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dhd1(t, D, x, prop=None, mRP=None):
    '''Compute partial derivatives of enthalpy w.r.t. t, p, or D at constant
    t, p, or D as a function of temperature, density, and composition'''
    def _rpfunc():
        return refprop.dhd1(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def fgcty(t, D, x, prop=None, mRP=None):
    '''Compute fugacity for each of the nc components of a mixture by
    numerical differentiation (using central differences) of the
    dimensionless residual Helmholtz energy'''
    def _rpfunc():
        return refprop.fgcty(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

#~ def fgcty2(t, D, x, prop=None, mRP=None):
    #~ '''Compute fugacity for each of the nc components of a mixture by
    #~ analytical differentiation of the dimensionless residual Helmholtz energy.'''
    #~ def _rpfunc():
        #~ return refprop.fgcty2(t, D, x)
    #~ return _rpfunc_handler(prop, mRP, _rpfunc)

def dbdt(t, x, prop=None, mRP=None):
    '''Compute the 2nd derivative of B (B is the second virial coefficient)
    with respect to T as a function of temperature and composition.'''
    def _rpfunc():
        return refprop.dbdt(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def virb(t, x, prop=None, mRP=None):
    '''Compute second virial coefficient as a function of temperature and
    composition.'''
    def _rpfunc():
        return refprop.virb(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def virc(t, x, prop=None, mRP=None):
    '''Compute the third virial coefficient as a function of temperature and
    composition.'''
    def _rpfunc():
        return refprop.virc(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def vird(t, x, prop=None, mRP=None):
    '''Compute the fourth virial coefficient as a function of temperature
    and composition.'''
    def _rpfunc():
        return refprop.vird(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def virba(t, x, prop=None, mRP=None):
    '''Compute second acoustic virial coefficient as a function of temperature
    and composition.'''
    def _rpfunc():
        return refprop.virba(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def virca(t, x, prop=None, mRP=None):
    '''Compute third acoustic virial coefficient as a function of temperature
    and composition.'''
    def _rpfunc():
        return refprop.virca(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def satt(t, x, kph=2, prop=None, mRP=None):
    '''Iterate for saturated liquid and vapor states given temperature and
    the composition of one phase'''
    def _rpfunc():
        return refprop.satt(t, x, kph)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def satp(p, x, kph=2, prop=None, mRP=None):
    '''Iterate for saturated liquid and vapor states given pressure and the
    composition of one phase.'''
    def _rpfunc():
        return refprop.satp(p, x, kph)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def satd(D, x, kph=2, prop=None, mRP=None):
    '''Iterate for temperature and pressure given a density along the
    saturation boundary and the composition.'''
    def _rpfunc():
        return refprop.satd(D, x, kph)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def sath(h, x, kph=2, prop=None, mRP=None):
    '''Iterate for temperature, pressure, and density given enthalpy along
    the saturation boundary and the composition.'''
    def _rpfunc():
        return refprop.sath(h, x, kph)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def sate(e, x, kph=2, prop=None, mRP=None):
    '''Iterate for temperature, pressure, and density given energy along the
    saturation boundary and the composition.'''
    def _rpfunc():
        return refprop.sate(e, x, kph)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def sats(s, x, kph=2, prop=None, mRP=None):
    '''Iterate for temperature, pressure, and density given entropy along
    the saturation boundary and the composition.'''
    def _rpfunc():
        return refprop.sats(s, x, kph)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def csatk(icomp, t, kph=2, prop=None, mRP=None):
    '''Compute the heat capacity along the saturation line as a function of
    temperature for a given component.'''
    def _rpfunc():
        return refprop.csatk(icomp, t, kph)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dptsatk(icomp, t, kph=2, prop=None, mRP=None):
    '''Compute the heat capacity and dP/dT along the saturation line as a
    function of temperature for a given component.'''
    def _rpfunc():
        return refprop.dptsatk(icomp, t, kph)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def cv2pk(icomp, t, D=0, prop=None, mRP=None):
    '''Compute the isochoric heat capacity in the two phase (liquid+vapor)
    region.'''
    def _rpfunc():
        return refprop.cv2pk(icomp, t, D)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def tprho(t, p, x, kph=2, kguess=0, D=0, prop=None, mRP=None):
    '''Iterate for density as a function of temperature, pressure, and
    composition for a specified phase.'''
    def _rpfunc():
        return refprop.tprho(t, p, x, kph, kguess, D)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def flsh(routine, var1, var2, x, kph=1, prop=None, mRP=None):
    '''Flash calculation given two independent variables and bulk
    composition.'''
    def _rpfunc():
        return refprop.flsh(routine, var1, var2, x, kph)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def flsh1(routine, var1, var2, x, kph=1, Dmin=0, Dmax=0, prop=None, mRP=None):
    '''Flash calculation given two independent variables and bulk
    composition.'''
    def _rpfunc():
        return refprop.flsh1(routine, var1, var2, x, kph, Dmin, Dmax)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def flsh2(routine, var1, var2, x, kq=1, ksat=0, tbub=0, tdew=0, pbub=0, pdew=0,
          Dlbub=0, Dvdew=0, xbub=None, xdew=None, prop=None, mRP=None):
    '''Flash calculation given two independent variables and bulk composition'''
    def _rpfunc():
        return refprop.flsh2(
            routine, var1, var2, x, kq, ksat, tbub, tdew, pbub, pdew, Dlbub,
            Dvdew, xbub, xdew)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def info(icomp=1, prop=None, mRP=None):
    '''Provides fluid constants for specified component.'''
    def _rpfunc():
        return refprop.info(icomp)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def xmass(x, prop=None, mRP=None):
    '''Converts composition on a mole fraction basis to mass fraction.'''
    def _rpfunc():
        return refprop.xmass(x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def xmole(xkg, prop=None, mRP=None):
    '''Converts composition on a mass fraction basis to mole fraction.'''
    def _rpfunc():
        return refprop.xmole(xkg)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def limitx(x, htype='EOS', t=0, D=0, p=0, prop=None, mRP=None):
    '''returns limits of a property model as a function of composition
    and/or checks input t, D, p against those limits.'''
    def _rpfunc():
        return refprop.limitx(x, htype, t, D, p)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def limitk(htype='EOS', icomp=1, t='tnbp', D=0, p=0, prop=None, mRP=None):
    '''Returns limits of a property model (read in from the .fld files) for
    a mixture component and/or checks input t, D, p against those limits.'''
    def _rpfunc():
        return refprop.limitk(htype, icomp, t, D, p)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def limits(x, htype = 'EOS', prop=None, mRP=None):
    '''Returns limits of a property model as a function of composition.'''
    def _rpfunc():
        return refprop.limits(x, htype)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def qmass(q, xliq, xvap, prop=None, mRP=None):
    '''converts quality and composition on a mole basis to a mass basis.'''
    def _rpfunc():
        return refprop.qmass(q, xliq, xvap)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def qmole(qkg, xlkg, xvkg, prop=None, mRP=None):
    '''Converts quality and composition on a mass basis to a molar basis.'''
    def _rpfunc():
        return refprop.qmole(qkg, xlkg, xvkg)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def wmol(x, prop=None, mRP=None):
    '''Molecular weight for a mixture of specified composition.'''
    def _rpfunc():
        return refprop.wmol(x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dielec(t, D, x, prop=None, mRP=None):
    '''Compute the dielectric constant as a function of temperature,
    density, and composition.'''
    def _rpfunc():
        return refprop.dielec(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def surft(t, x, prop=None, mRP=None):
    '''Compute surface tension.'''
    def _rpfunc():
        return refprop.surft(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def surten(t, Dliq, Dvap, xliq, xvap, prop=None, mRP=None):
    '''Compute surface tension.'''
    def _rpfunc():
        return refprop.surten(t, Dliq, Dvap, xliq, xvap)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def meltt(t, x, prop=None, mRP=None):
    '''Compute the melting line pressure as a function of temperature and
    composition.'''
    def _rpfunc():
        return refprop.meltt(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def meltp(p, x, prop=None, mRP=None):
    '''Compute the melting line temperature as a function of pressure and
    composition.'''
    def _rpfunc():
        return refprop.meltp(p, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def sublt(t, x, prop=None, mRP=None):
    '''Compute the sublimation line pressure as a function of temperature
    and composition.'''
    def _rpfunc():
        return refprop.sublt(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def sublp(p, x, prop=None, mRP=None):
    '''Compute the sublimation line temperature as a function of pressure
    and composition.'''
    def _rpfunc():
        return refprop.sublp(p, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def trnprp(t, D, x, prop=None, mRP=None):
    '''Compute the transport properties of thermal conductivity and
    viscosity as functions of temperature, density, and composition.'''
    def _rpfunc():
        return refprop.trnprp(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def getktv(icomp, jcomp, prop=None, mRP=None):
    '''Retrieve mixture model and parameter info for a specified binary.'''
    def _rpfunc():
        return refprop.getktv(icomp, jcomp)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def setktv(icomp, jcomp, hmodij, fij=([0] * refprop._nmxpar),
           hfmix='hmx.bnc'):
    '''Set mixture model and/or parameters.'''
    _checksetupblock('setktv')
    return refprop.setktv(icomp, jcomp, hmodij, fij, hfmix)

def setaga():
    '''Set up working arrays for use with AGA8 equation of state.'''
    _checksetupblock('setaga')
    return refprop.setaga()

def unsetaga():
    '''Load original values into arrays changed in the call to SETAGA.'''
    _checksetupblock('unsetaga')
    return refprop.unsetaga()

def preos(ixflag):
    '''Turn on or off the use of the PR cubic equation.'''
    _checksetupblock('preos')
    return refprop.preos(ixflag)

def getfij(hmodij, prop=None, mRP=None):
    '''Retrieve parameter info for a specified mixing rule.'''
    def _rpfunc():
        return refprop.getfij(hmodij)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def b12(t, x, prop=None, mRP=None):
    '''Compute b12 as a function of temperature and composition.'''
    def _rpfunc():
        return refprop.b12(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def excess(t, p, x, kph=0, prop=None, mRP=None):
    '''Compute excess properties as a function of temperature, pressure, and
    composition.'''
    def _rpfunc():
        return refprop.excess(t, p, x, kph)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def getmod(icomp, htype, prop=None, mRP=None):
    '''Retrieve citation information for the property models used'''
    def _rpfunc():
        return refprop.getmod(icomp, htype)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def cstar(t, p, v, x, prop=None, mRP=None):
    '''Calculate the critical flow factor, C*, for nozzle flow of a gas
    (subroutine was originally named CCRIT)'''
    def _rpfunc():
        return refprop.cstar(t, p, v, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def phiderv(icomp, jcomp, t, D, x, prop=None, mRP=None):
    '''Calculate various derivatives needed for VLE determination'''
    def _rpfunc():
        return refprop.phiderv(icomp, jcomp, t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def chempot(t, D, x, prop=None, mRP=None):
    '''Compute the chemical potentials for each of the nc components of a
    mixture.'''
    def _rpfunc():
        return refprop.chempot(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def fugcof(t, D, x, prop=None, mRP=None):
    '''Compute the fugacity coefficient for each of the nc components of a
    mixture.'''
    def _rpfunc():
        return refprop.fugcof(t, D, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dcdt(t, x, prop=None, mRP=None):
    '''Compute the 1st derivative of C (C is the third virial coefficient) with
    respect to T as a function of temperature and composition.'''
    def _rpfunc():
        return refprop.dcdt(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def dcdt2(t, x, prop=None, mRP=None):
    '''Compute the 2nd derivative of C (C is the third virial coefficient) with
    respect to T as a function of temperature and composition.'''
    def _rpfunc():
        return refprop.dcdt2(t, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def fpv(t, D, p, x, prop=None, mRP=None):
    '''Compute the supercompressibility factor, Fpv.'''
    def _rpfunc():
        return refprop.fpv(t, D, p, x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

def rmix2(x, prop=None, mRP=None):
    '''Return the gas "constant" as a combination of the gas constants for
    the pure fluids.'''
    def _rpfunc():
        return refprop.rmix2(x)
    return _rpfunc_handler(prop, mRP, _rpfunc)

#compilations

#create function (to demonstrate pipe)
def _mRPpipe(fluid, mRP):
    resetup(fluid, mRP=mRP)
    for each in range(100):
        #put value into pipe
        mRP['cpipe'].send(
            press(303 + each, 58, [0.4, 0.6])['p'])
    mRP['cpipe'].close()

if __name__ == '__main__':
    #add module file path to python sys path
    import multiRP as _filename
    _filename = (os.path.dirname(_filename.__file__))
    sys.path.append(_filename)

    #test multiRP without multiprocessing
    import rptest
    rptest.settest('multiRP')

    #initiate multirefprop, this will create global 'mRP'
    multirefprop()

    #create setup details
    H2O = setup('def', 'WATER')
    H2O_NH3 = setup('def', 'AMMONIA', 'WATER')
    CH4 = setup('def', 'METHANE')
    CH4_C2H6 = setup('def', 'METHANE', 'ETHANE')

    #various refprop functions named
    peen = mRP['process'](target=critp, args=([1], H2O, mRP))
    ptwee = mRP['process'](target=critp, args=([0.4, 0.6],),
                           kwargs={'prop':H2O_NH3, 'mRP':mRP})#alternative input
    pdrie = mRP['process'](target=critp, args=([1], CH4, mRP))
    pvier = mRP['process'](target=critp, args=([0.35, 0.65], CH4_C2H6, mRP))
    qeen = mRP['process'](target=critp, args=([1], H2O, mRP))
    qtwee = mRP['process'](target=critp, args=([0.4, 0.6], H2O_NH3, mRP))
    qdrie = mRP['process'](target=critp, args=([1], CH4, mRP))
    qvier = mRP['process'](target=critp, args=([0.35, 0.65], CH4_C2H6, mRP))
    ween = mRP['process'](target=critp, args=([1], H2O, mRP))
    wtwee = mRP['process'](target=critp, args=([0.4, 0.6], H2O_NH3, mRP))
    wdrie = mRP['process'](target=critp, args=([1], CH4, mRP))
    wvier = mRP['process'](target=critp, args=([0.35, 0.65], CH4_C2H6, mRP))
    reen = mRP['process'](target=critp, args=([1], H2O, mRP))
    rtwee = mRP['process'](target=critp, args=([0.4, 0.6], H2O_NH3, mRP))
    rdrie = mRP['process'](target=critp, args=([1], CH4, mRP))
    rvier = mRP['process'](target=critp, args=([0.35, 0.65],
                                               CH4_C2H6, mRP))
    sfun = mRP['process'](target=_mRPpipe, args=(H2O_NH3, mRP))#pipe

    #list refprop functions
    processlist = [peen, ptwee, pdrie, pvier, qeen, qtwee, qdrie, qvier,
                   ween, wtwee, wdrie, wvier, reen, rtwee, rdrie, rvier,
                   sfun]

    #now run the functions
    run_mRP(processlist)

    #print the results
    for each in processlist:
        if each.name in mRP['result']:
            print(mRP['result'][each.name]) #only returns the last refprop function results unless tweaked at function call

    #loop untill pipe is empty
    while mRP['ppipe'].poll():#parentpipe
        #print value from pipe
        print(mRP['ppipe'].recv())
