#-------------------------------------------------------------------------------
#Name:                  equipment
#Purpose:               Calculate equipment fluid output properties
#                       using refprop module
#
#Author:                Thelen, B.J.
#                       thelen_ben@yahoo.com
#-------------------------------------------------------------------------------

u'''This module usage the refprop module to calculate fluid output
properties from various mechanical equipment.

Refer to refprop module for details on calling fluids and fluid properties.

Calculated fluid output properties will be added into the refprop 'prop'
library. Equipment characteristics will be added into a 'equip' library.

Units
----------------------------------------------------------------------------
temperature                             K
pressure, fugacity                      kPa
density                                 mol/L
composition                             mole fraction
quality                                 mole basis (moles vapor/total moles)
enthalpy, internal energy               J/mol
Gibbs, Helmholtz free energy            J/mol
entropy, heat capacity                  J/(mol.K)
speed of sound                          m/s
Joule-Thompson coefficient              K/kPa
d(p)/d(rho)                             kPa.L/mol
d2(p)/d(rho)2                           kPa.(L/mol)^2
viscosity                               microPa.s (10^-6 Pa.s)
thermal conductivity                    W/(m.K)
dipole moment                           debye
surface Tension                         N/m
composition flow                        mole/s
----------------------------------------------------------------------------'''
from __future__ import division
import os, sys, scipy, math, decimal, copy
import multiRP
import multiprocessing as mp

#classes
class EquipmentError(Exception):
    u'General EquipmentError for python module'
    pass

class EquipmentinputError(EquipmentError):
    u'Equipment input Error'
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class EquipmentinfeasibleError(EquipmentError):
    u'Equipment input Error'
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class SetMultiprocessing(object):
    u'Return _setmultiprocessing status (on / off)'
    def __repr__(self):
        if not u'_setmultiprocessing' in globals(): SetMultiprocessing.on()
        return _setmultiprocessing
    @staticmethod
    def on():
        u'Sets SetMultiprocessing on, allows equipment multiprocessing'
        global _setmultiprocessing
        _setmultiprocessing = u'on'
        return _setmultiprocessing
    @staticmethod
    def off():
        u'Sets SetMultiprocessing off, single processing for equipment module'
        global _setmultiprocessing
        _setmultiprocessing = u'off'
        return _setmultiprocessing


#suporting functions
def _inputerrorcheck(deflocals, equipment):
    u'''check valid def input or else raise error'''
    for key in deflocals.keys():
        #check input prop
        if key[:4] == u'prop':
            if deflocals[key] != None:
                if not deflocals[key].__class__ == dict:
                    raise EquipmentinputError(u"expect dict input for 'prop' " +
                                               u"instead of " +
                                               unicode(deflocals[u'prop'].__class__))

            #check input prop['t']
            if [u'subcool', u'exchanger'].__contains__(equipment):
                if deflocals[key] != None:
                    if not u't' in deflocals[key]:
                        raise EquipmentinputError(u'prop dict key "t" ' +
                                                   u'required as input')

            #check input prop['h']
            if [u'separator', u'reg_valve', u'turbine', u'flowmerge',
                u'turbine_separator',
                u'turbine_bleed_sep'].__contains__(equipment):
                if deflocals[key] != None:
                    if not u'h' in deflocals[key]:
                        raise EquipmentinputError(u'prop dict key "h" ' +
                                                   u'required as input')

            #check input prop['x']
            if [u'subcool', u'pump', u'separator', u'reg_valve', u'turbine',
                u'flowmerge', u'turbine_separator', u'compressor',
                u'turbine_bleed_sep', u'exchanger'].__contains__(equipment):
                if deflocals[key] != None:
                    if not u'x' in deflocals[key]:
                        raise EquipmentinputError(u'prop dict key "x" ' +
                                                   u'required as input')

            #check input prop['Q']
            if [u'pump', u'separator', u'turbine', u'flowmerge',
            u'turbine_separator', u'compressor', u'turbine_bleed_sep',
            u'exchanger'].__contains__(equipment):
                if deflocals[key] != None:
                    if not u'Q' in deflocals[key]:
                        raise EquipmentinputError(u'prop dict key "Q" ' +
                                                   u'required as input')
                    elif deflocals[key][u'Q'] < 0:
                        raise EquipmentinputError(u'"Q" value to be > 0')

            #check input prop['s']
            if [u'pump', u'turbine', u'compressor'].__contains__(equipment):
                if deflocals[key] != None:
                    if not u's' in deflocals[key]:
                        raise EquipmentinputError(u'prop dict key "s" ' +
                                                   u'required as input')

            #check input prop['p']
            if [u'pump', u'separator',u'flowmerge', u'turbine_separator',
                u'compressor', u'exchanger',
                u'turbine_bleed_sep'].__contains__(equipment):
                if deflocals[key] != None:
                    if not u'p' in deflocals[key]:
                        raise EquipmentinputError(u'prop dict key "p" ' +
                                                   u'required as input')

            #check equal prop_1 and prop_2['p']
            if [u'flowmerge'].__contains__(equipment):
                if u'prop_2' in deflocals.keys() \
                and deflocals[u'prop_1'][u'p'] != deflocals[u'prop_2'][u'p']:
                    raise EquipmentinputError(u'pressure value to be equal ' +
                                               u'for prop_1 and prop_2')

            #check equal prop_1 and prop_2['hfld']
            if [u'flowmerge'].__contains__(equipment):
                if u'prop_2' in deflocals.keys() \
                and deflocals[u'prop_1'][u'hfld'] != deflocals[u'prop_2'][u'hfld']:
                    raise EquipmentinputError(u'hfld value to be equal for ' +
                                               u'prop_1 and prop_2')

        #check input dt
        if ([u'subcool'].__contains__(equipment) \
        or [u'exchanger'].__contains__(equipment)) \
        and key == u'dt':
            if not deflocals[key].__class__ == float \
            and not deflocals[key].__class__ == int:
                raise EquipmentinputError(u"expect float/int input for 'dt' " +
                                           u'instead of ' +
                                           unicode(deflocals[key].__class__))
            if deflocals[u'dt'] < 0:
                raise EquipmentinputError(u"'dt' value shall be positive")

        #check input eff
        if [u'pump', u'turbine', u'turbine_separator', u'turbine_bleed_sep',
            u'compressor'].__contains__(equipment) \
            and (key == u'eff' \
                 or key == u'e_eff'):
            if not deflocals[key].__class__ == float \
            and not 0 <= key <= 1:
                raise EquipmentinputError(u"expect float value between 0 and 1 "
                                           u"for '" + key + u"' instead of " +
                                           unicode(deflocals[key]))

        #check input Q_ratio
        if [u'turbine_bleed_sep'].__contains__(equipment) \
           and (key == u'Q_ratio' or key == u'Q_ratiobleed'):
            if not deflocals[key].__class__ == float \
            and not 0 < key < 1:
                raise EquipmentinputError(u"expect float value between 0 and 1 "
                                           u"for '" + key + u"' instead of " +
                                           unicode(deflocals[key]))

        #check input dp
        if [u'pump', u'separator', u'turbine_separator', u'compressor', u'exchanger',
            u'turbine_bleed_sep'].__contains__(equipment) and key == u'dp':
            if not deflocals[key].__class__ == float \
            and not deflocals[key].__class__ == int:
                raise EquipmentinputError(u"expect float/int input for 'dp' " +
                                           u'instead of ' +
                                           unicode(deflocals[key].__class__))
            if deflocals[u'dp'] < 0:
                raise EquipmentinputError(u"'dp' value shall be positive")

        #check input dp_contra
        if [u'exchanger'].__contains__(equipment) and key == u'dp_contra':
            if not deflocals[key].__class__ == float \
            and not deflocals[key].__class__ == int:
                raise EquipmentinputError(u"expect float/int input for '" +
                                           u"'dp_contra' instead of " +
                                           unicode(deflocals[key].__class__))
            if deflocals[u'dp_contra'] < 0:
                raise EquipmentinputError(u"'dp_contra' value shall be positive")

        #check input p_out
        if [u'reg_valve', u'turbine', u'turbine_bleed_sep',
             u'turbine_separator'].__contains__(equipment) \
             and (key == u'p_out' or key == u'p_outbleed'):
            if not deflocals[key].__class__ == float \
            and not deflocals[key].__class__ == int:
                raise EquipmentinputError(u"expect float/int input for '" +
                                           key + u"' instead of " +
                                           unicode(deflocals[key].__class__))


def _prop(name, value, unit, valueSI=None, unitSI=None):
    u"""print out fluid property in organised format

    input:
        name--fluid characteristic type e.g. pressure
        value--characteristic value
        unit--characteristic unit e.g. kPa
        ValueSI--characteristic value in SI units
        unitSI--characteristic SI unit e.g. bar
    output:
        fluid object"""
    if (value == scipy.inf or value == -scipy.inf) and valueSI == None:
        #format print out
        fldob = unicode(u'{0:22}' + u'{1:>13,.3f}' + u'{2:18}\n').format(name, value,
                                                                  u'')
    if valueSI == None:
        #format print out
        fldob = unicode(u'{0:22}' + u'{1:>13,.3f}' + u'{2:18}\n').format(name, value,
                                                                  unit)
    elif value == scipy.inf or value == -scipy.inf:
        #format print out
        fldob = unicode(u'{0:22}' + u'{1:>13,.3f}' + u'{2:19}' + u'{3:>13,.3f}' +
                    u'{4}\n').format(name, value, u'', valueSI, u'')
    else:
        #format print out
        fldob = unicode(u'{0:22}' + u'{1:>13,.3f}' + u'{2:18}' + u'{3:>13,.3f}' +
                    u'{4}\n').format(name, value, unit, valueSI, unitSI)
    return fldob


def _fluidprop(prop_in):
    u"""print all fluid properties in organised format

    input:
        prop_in--fluid characteristics to be plotted
    output:
        prpstr--fluid property in organized string format"""
    prpstr = u''
    if u'hfld' in prop_in.keys():
        multiRP.resetup(prop_in)
        if u'x' in prop_in.keys():
            x = []
            for each in xrange(len(prop_in[u'x'])):
                x.insert(each, float(prop_in[u'x'][each]))
        if u'xliq' in prop_in.keys():
            xliq = []
            for each in xrange(len(prop_in[u'xliq'])):
                xliq.insert(each, float(prop_in[u'xliq'][each]))
        if u'xvap' in prop_in.keys():
            xvap = []
            for each in xrange(len(prop_in[u'xvap'])):
                xvap.insert(each, float(prop_in[u'xvap'][each]))
        #determine mole weigth composition
        #assigned [0] by echanger to correct for wmix = 0
        if x == [0]:
            wmix = 0
        else:
            wmix = multiRP.xmass(x)[u'wmix']
        if u'q' in prop_in.keys():
            if 0 < prop_in[u'q'] < 1:
                if u'xliq' in prop_in.keys():
                    wmixliq = multiRP.xmass(xliq)[u'wmix']
                if u'xvap' in prop_in.keys():
                    wmixvap = multiRP.xmass(xvap)[u'wmix']
        #fluid components display
        #correction on exchanger assigned phantom fluid
        if u'hfld' in prop_in.keys():
            for each in xrange(len(prop_in[u'hfld'])):
                #obtain mole weigth individual comp.
                wmm = multiRP.info(each + 1)[u'wmm']
                prpstr += _prop(unicode(prop_in[u'hfld'][each]),
                                x[each] * 100,
                                u'% mol(i)/mol',
                                x[each] / wmix * wmm * 100,
                                u'% g(i)/g')
        multiRP.resetup(prop_in)
        #liquid components display
        if u'q' in prop_in.keys():
            if 0 < prop_in[u'q'] < 1:
                for each in xrange(len(prop_in[u'hfld'])):
                    #obtain mole weigth individual comp.
                    wmm = multiRP.info(each + 1)[u'wmm']
                    prpstr += _prop(unicode(prop_in[u'hfld'][each]) + u' (liq)',
                                    xliq[each] * 100,
                                    u'% mol(i)/mol',
                                    xliq[each] / wmixliq * wmm * 100,
                                    u'% g(i)/g')
        multiRP.resetup(prop_in)
        #vapor components display
        if u'q' in prop_in.keys():
            if 0 < prop_in[u'q'] < 1:
                for each in xrange(len(prop_in[u'hfld'])):
                    #obtain mole weigth individual comp.
                    wmm = multiRP.info(each + 1)[u'wmm']
                    prpstr += _prop(unicode(prop_in[u'hfld'][each]) + u' (vap)',
                                    xvap[each] * 100,
                                    u'% mol(i)/mol',
                                    xvap[each] / wmixvap * wmm * 100,
                                    u'% g(i)/g')
        multiRP.resetup(prop_in)
        #Flow display
        if u'Q' in prop_in.keys():
            prpstr += _prop(u'Flow: ',
                            prop_in[u'Q'],
                            u'  mol/s',
                            prop_in[u'Q'] * wmix * 3600 / 1000,
                            u'  kg/h')
        #mol weight display
        #correction for exchanger phantom fluid
        if not wmix == 0:
            prpstr += _prop(u'molecular weight: ',
                            wmix,
                            u'  g/mol')
        #pressure display
        if u'p' in prop_in.keys():
            prpstr += _prop(u'pressure: ',
                            prop_in[u'p'],
                            u'  kPa(a)',
                            prop_in[u'p'] / 100,
                            u'  bar(a)')
        #temp display
        if u't' in prop_in.keys():
            prpstr += _prop(u'temperature: ',
                            prop_in[u't'],
                            u'  Kelvin',
                            prop_in[u't'] - 273.15,
                            u'  Celsius')
        #Density display
        if u'D' in prop_in.keys():
            prpstr += _prop(u'Density: ',
                            prop_in[u'D'],
                            u'  mol/L',
                            prop_in[u'D'] * wmix,
                            u'  kg/m3')
        #liq Density display
        if u'q' in prop_in.keys():
            if 0 < prop_in[u'q'] < 1:
                if u'Dliq' in prop_in.keys():
                    prpstr += _prop(u'Density (liq): ',
                                    prop_in[u'Dliq'],
                                    u'  mol/L',
                                    prop_in[u'Dliq'] * wmixliq,
                                    u'  kg/m3')
        #vap density display
        if u'q' in prop_in.keys():
            if 0 < prop_in[u'q'] < 1:
                if u'Dvap' in prop_in.keys():
                    prpstr += _prop(u'Density (vap): ',
                                    prop_in[u'Dvap'],
                                    u'  mol/L',
                                    prop_in[u'Dvap'] * wmixvap,
                                    u'  kg/m3')
        #enthalpy display
        if u'h' in prop_in.keys() and prop_in[u'h'] != 0:
            prpstr += _prop(u'Enthalpy: ',
                            prop_in[u'h'],
                            u'  J/mol',
                            prop_in[u'h'] / wmix,
                            u'  kJ/kg')
        #entropy display
        if u's' in prop_in.keys() and prop_in[u'h'] != 0:
            prpstr += _prop(u'Entropy: ',
                            prop_in[u's'],
                            u'  J/(molK)',
                            prop_in[u's'] / wmix,
                            u'  kJ/kgK')
        #phase display
        try:
            phase = multiRP.getphase(prop_in)
            #quality display
            if phase == u"2 phase":
                prpstr += _prop(u'vapor quality: ', prop_in[u'q'], u'  mol(v)/mol',
                                prop_in[u'q'] * wmixvap / wmix, u'  g(v)/g')
            #liquid (subcool) display
            elif phase == u"liquid":
                t = multiRP.satp(prop_in[u'p'], x, 1)[u't']
                prpstr += _prop(u'subcool: ', t - prop_in[u't'], u'  Kelvin',
                                t - prop_in[u't'], u'  Celsius')
            #saturated liquid display
            elif phase == u"saturated liquid":
                prpstr += u'saturated liquid\n'
            #saturated vapor display
            elif phase == u"saturated vapor":
                prpstr += u'saturated vapor\n'
            #superheated vapor display
            elif phase == u"vapor":
                t = multiRP.satp(prop_in[u'p'], prop_in[u'x'], 2)[u't']
                prpstr += _prop(u'superheated vapor: ', prop_in[u't'] - t, u'  Kelvin',
                                prop_in[u't'] - t, u'  Celsius')
            #superheated gas display
            elif phase == u"gas":
                t = multiRP.satp(prop_in[u'p'], prop_in[u'x'], 2)[u't']
                prpstr += _prop(u'superheated gas: ', prop_in[u't'] - t, u'  Kelvin',
                                prop_in[u't'] - t, u'  Celsius')
            #compressible liquid display
            elif phase == u"compressible liquid":
                pcrit = multiRP.critp(prop_in[u'x'])[u'pcrit']
                prpstr += _prop(u'compres. lqd: ', prop_in[u'p'] - pcrit,
                                u'  kPa(a)', (prop_in[u'p'] - pcrit) / 100,
                                u'  bar(a)')
            #supercritical fluid display
            elif phase == u"Supercritical fluid":
                prop = multiRP.critp(prop_in[u'x'])
                pcrit = prop[u'pcrit']
                tcrit = prop[u'tcrit']
                prpstr += _prop(u'Supercritical: ', prop_in[u'p'] - pcrit,
                                u'  kPa(a)', (prop_in[u'p'] - pcrit) / 100,
                                u'  bar(a)')
                prpstr += _prop(u'', prop_in[u't'] - tcrit, u'  Kelvin',
                                prop_in[u't'] - tcrit, u'  Celsius')
        except:
            pass
    #return fluid object
    return prpstr


def _exchgraph(prop_in, prop_out, contra_prop_in, contra_prop_out, mRP=None):
    u"""return a string representing the exchanger in /out fluids temp graph"""
    try:
        lines = 30
        position = 80
        c, h = [], []
        #define cold and hot fluids
        if prop_in[u't'] < contra_prop_in[u't']:
            cold_in = prop_in
            cold_out = prop_out
            hot_in = contra_prop_in
            hot_out = contra_prop_out
        else:
            cold_in = contra_prop_in
            cold_out = contra_prop_out
            hot_in = prop_in
            hot_out = prop_out
        #calculate temperature at 80 positions (line length = 80)
        #determin h and p step size
        h_cold_step = (cold_out[u'h'] - cold_in[u'h'])/(position-1)
        p_cold_step = (cold_out[u'p'] - cold_in[u'p'])/(position-1)
        h_hot_step = (hot_in[u'h'] - hot_out[u'h'])/(position-1)
        p_hot_step = (hot_in[u'p'] - hot_out[u'p'])/(position-1)

        #check if multiprocessing is allowed
        if SetMultiprocessing().__repr__() == u'on':
            ##start with multiprocessing##
            if mRP == None:
                #create if not exist
                mRP = multiRP.multirefprop()#initiate multirefprop

            ##create children
            for pos in xrange(position):
                c.append(mRP[u'process'](target=_coldhot,
                                        args=(cold_in[u'p'] + p_cold_step * pos,
                                              cold_in[u'h'] + h_cold_step * pos,
                                              cold_in[u'x']),
                                        kwargs={u'prop':cold_in, u'mRP':mRP}))
                h.append(mRP[u'process'](target=_coldhot,
                                        args=(hot_out[u'p'] + p_hot_step * pos,
                                              hot_out[u'h'] + h_hot_step * pos,
                                              hot_in[u'x']),
                                        kwargs={u'prop':hot_in, u'mRP':mRP}))
            fldlist = c + h

            #run the multiprocessing list
            multiRP.run_mRP(fldlist)
            #obtain results from multiprocessing list
            templist = []
            for each in fldlist:
                if each.name in mRP[u'result']:
                    if mRP[u'result'][each.name].__class__ == dict \
                    and u't' in mRP[u'result'][each.name]:
                        templist.append(mRP[u'result'][each.name][u't'])
                    else:
                        templist.append(mRP[u'result'][each.name])
            ##end with multiprocessing##

        #check for single processing
        elif SetMultiprocessing().__repr__() == u'off':
            #loop to get each position
            for pos in xrange(position):
                c.append(_coldhot(cold_in[u'p'] + p_cold_step * pos,
                                  cold_in[u'h'] + h_cold_step * pos,
                                  cold_in[u'x'],
                                  prop=cold_in))
                h.append(_coldhot(hot_out[u'p'] + p_hot_step * pos,
                                  hot_out[u'h'] + h_hot_step * pos,
                                  hot_in[u'x'],
                                  prop=hot_in))
            #create complete temp. list for each position
            templist = [each if each == None else each[u't'] for each in (c + h)]

        #round of due to graph display hickups
        for each in xrange(len(templist)):
            if templist[each] != None:
                templist[each] = round(templist[each], 10)

        #split templist (hot and cold)
        c = templist[:int(len(templist)/2)]
        h = templist[int(len(templist)/2):]

        #determine temperature range and split evenly for 30 lines
        dt_step = (hot_in[u't'] - cold_in[u't']) / (lines - 1)
        #create lines
        graph = u''
        for line in xrange(lines):
            #create positions within line
            for pos in xrange(position):
                #if both h and c then plot X
                if c[pos] != None and h[pos] != None \
                and hot_in[u't'] - dt_step * (line + 1) < \
                c[pos] <= \
                hot_in[u't'] - dt_step * line \
                and hot_in[u't'] - dt_step * (line + 1) < \
                h[pos] <= \
                hot_in[u't'] - dt_step * line:
                    graph += u'X'
                #if c only then plot C
                elif c[pos] != None \
                and hot_in[u't'] - dt_step * (line + 1) < \
                c[pos] <= \
                hot_in[u't'] - dt_step * line:
                    graph += u'C'
                #if h only then plot H
                elif h[pos] != None \
                and hot_in[u't'] - dt_step * (line + 1) < \
                h[pos] <= \
                hot_in[u't'] - dt_step * line:
                    graph += u'H'
                #else plot space
                else: graph += u' '
            #next line
            graph += u'\n'
        return graph
    except:
        return u'error in generating graphics at def _exchgraph\n'

def _coldhot(p, h, x, prop=None, mRP=None):
    u'''cold / hot function for _exchgraph'''
    try:
        return multiRP.flsh(u'ph', p, h, x, prop=prop, mRP=mRP)
    except multiRP.RefpropError:
        #check for multiprocessing
        if SetMultiprocessing().__repr__() == u'on':
            multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = None
        #check for single processing
        elif SetMultiprocessing().__repr__() == u'off':
            return None



#primary functions
def pump(prop_in, dp, eff=1.0, e_eff=1.0, name=False, mRP=None):
    u"""Calculate fluid properties at pump discharge.
    Equipment specification output through Pump

    input:
        prop_in--fluid properties library contain keys:
            'Q'--mole flow (mole / sec)
            's'--entropy
            'p'--pressure
            'x'--mole fraction
            'h'--enthalpy (optional)
        dp--pump head
        eff--pump efficiency (100% std)
        e_eff--e-motor efficiency (100% std)
    output:
        prop--fluid properties
        equipspec--equipment details"""
    #global fluid, Pump
    _inputerrorcheck(locals(), u'pump')
    multiRP.resetup(prop_in, mRP=mRP)
    p = prop_in[u'p'] + dp
    s = prop_in[u's']
    x = prop_in[u'x']
    flshcalc = False
    if u'h' in prop_in: h = prop_in[u'h']
    else:
        #inhibit error reporting
        if unicode(multiRP.SetErrorDebug()) == u'on':
            sed = multiRP.SetErrorDebug.on
            multiRP.SetErrorDebug.off()
        elif unicode(multiRP.SetErrorDebug()) == u'off':
            sed = multiRP.SetErrorDebug.off
        try: #use flsh1 calculation for speed increase
            h = multiRP.psliq(p, s, x)[u'h']
            prop.update(sed())
        except multiRP.RefpropError:
            sed()
            h = multiRP.flsh(u'ps', p, s, x)[u'h']
            flshcalc = True
    Q = prop_in[u'Q']

    #calculate enthalpy at 100% efficiency
    #inhibit error reporting
    if unicode(multiRP.SetErrorDebug()) == u'on':
        sed = multiRP.SetErrorDebug.on
        multiRP.SetErrorDebug.off()
    elif unicode(multiRP.SetErrorDebug()) == u'off':
        sed = multiRP.SetErrorDebug.off
    try: #use flsh1 calculation for speed increase
        if flshcalc: raise multiRP.RefpropError
        prop = multiRP.psliq(p, s, x)
        prop.update(sed())
    except multiRP.RefpropError:
        sed()
        prop = multiRP.flsh(u'ps', p, s, x)
        flshcalc = True

    #correct h for efficiency losses
    h += (prop[u'h'] - h) / eff

    #calculate properties with corrected h and discharge p
    #inhibit error reporting
    if unicode(multiRP.SetErrorDebug()) == u'on':
        sed = multiRP.SetErrorDebug.on
        multiRP.SetErrorDebug.off()
    elif unicode(multiRP.SetErrorDebug()) == u'off':
        sed = multiRP.SetErrorDebug.off
    try: #use flsh1 calculation for speed increase
        if flshcalc: raise multiRP.RefpropError
        prop = multiRP.phliq(p, h, x)
        prop.update(sed())
    except multiRP.RefpropError:
        sed()
        prop = multiRP.flsh(u'ph', p, h, x)

    #restore Q value
    prop[u'Q'] = Q

    #power consumption
    pwr = ((prop[u'h'] - prop_in[u'h']) * Q) / e_eff

    equipspec ={u'type': u'PUMP', u'name': name, u'dp': dp, u'eff': eff,
                u'fluidin': prop_in, u'fluidout': prop, u'pwr': -pwr,
                u'e_eff':e_eff}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop


def exchanger(prop_in, dp, dt, prop_contra_in={}, dp_contra=0,
              prop_contra_out={}, name=False, mRP=None):
    u"""Calculate fluid out properties of exchanger.
    Equipment specification output through Exchanger

    input:
        prop_in--fluid properties library contain keys:
            'Q'--mole flow (mole / sec)
            't'--temperature
            'p'--pressure
            'x'--mole fraction
            'h'--enthalpy (optional)
            'q'--vapor component (optional)
        dp--delta pressure
        dt--min. dif. temperature between counterflow in/outlet
        prop_contra_in--fluid properties library contain keys (optional):
            'Q'--mole flow (mole / sec) (not req. to initiate calc.)
            't'--temperature (not req. to initiate calc.)
            'p'--pressure (not req. to initiate calc.)
            'x'--mole fraction (not req. to initiate calc.)
            'h'--enthalpy (optional)
            'q'--vapor component (optional)
        dp_contra--delta pressure contra flow (optional)
        prop_contra_out--fluid properties library contain keys:
            't'--temperature (not req. to initiate calc.)
            'h'--enthalpy (optional)
        cmplx--turn on / off complex calculations
            False--simple calculation assuming set value for dt only
            True--calculation also performs correction to prevent crossing of
                incomming and outcomming temperature lines.
    output:
        prop--fluid properties of calculated stream out
        equipspec--equipment details"""

    #create prop_contra_out equal to prop_contra_in
    #set prop_contra_out['t'] to 0 if not specified 'allows to initiate calc.'
    if (not u't' in locals()[u'prop_contra_out']
         or not u'Q' in locals()[u'prop_contra_out']
         or prop_contra_out[u'Q'] <= 0):
        prop_contra_out[u't'] = 0
    if not u'Q' in locals()[u'prop_contra_in']:
        prop_contra_out[u'Q'] = 0
        prop_contra_in[u'Q'] = 0
    else: prop_contra_out[u'Q'] = prop_contra_in[u'Q']
    if not u'p' in locals()[u'prop_contra_in']:
        prop_contra_out[u'p'] = 0
        prop_contra_in[u'p'] = 0
    else: prop_contra_out[u'p'] = prop_contra_in[u'p'] - dp_contra
    if not u'x' in locals()[u'prop_contra_in']:
        prop_contra_out[u'x'] = [0]
        prop_contra_in[u'x'] = [0]
    else: prop_contra_out[u'x'] = prop_contra_in[u'x']
    if not u't' in locals()[u'prop_contra_in']:
        prop_contra_in[u't'] = 0
    if not u'h' in locals()[u'prop_contra_in']:
        prop_contra_in[u'h'] = 0
    if not u'h' in locals()[u'prop_contra_out']:
        prop_contra_out[u'h'] = 0

    _inputerrorcheck(locals(), u'exchanger')

    #define variables with values
    Q_in = prop_in[u'Q']
    p_in = prop_in[u'p']
    t_in = prop_in[u't']
    x_in = prop_in[u'x']
    if not u'h' in prop_in:
        multiRP.resetup(prop_in, mRP=mRP)
        ########################################################################
        ############ actual h should be present for pure fluids ################
        ########################################################################
        h_in = multiRP.flsh(u'tp', t_in, p_in, x_in)[u'h']
    else: h_in = prop_in[u'h']
    p_out = p_in - dp
    if p_out < 0: p_out = 0

    #define 0 contra flow 'for initiating of calc'
    if prop_contra_out[u't'] == 0 or abs(t_in - prop_contra_in[u't']) < dt:
        multiRP.resetup(prop_in, mRP=mRP)
        prop = multiRP.flsh(u'ph', p_out, h_in, x_in)

    #define Dh based on calculations
    else:
        #define cold and hot prop stream
        if prop_in[u't'] <= prop_contra_in[u't']:
            prop_cold_in = copy.copy(prop_in)
            prop_hot_in = copy.copy(prop_contra_in)
            dp_cold = dp
            dp_hot = dp_contra
        elif prop_in[u't'] > prop_contra_in[u't']:
            prop_hot_in = copy.copy(prop_in)
            prop_cold_in = copy.copy(prop_contra_in)
            dp_cold = dp_contra
            dp_hot = dp
        Q_cold_in = prop_cold_in[u'Q']
        t_cold_in = prop_cold_in[u't']
        p_cold_in = prop_cold_in[u'p']
        x_cold_in = prop_cold_in[u'x']
        if not u'q' in prop_cold_in: q_cold_in = -1
        else: q_cold_in = prop_cold_in[u'q']
        if not u'h' in prop_cold_in:
            multiRP.resetup(prop_cold_in, mRP=mRP)
        else: h_cold_in = prop_cold_in[u'h']
        Q_hot_in = prop_hot_in[u'Q']
        t_hot_in = prop_hot_in[u't']
        p_hot_in = prop_hot_in[u'p']
        x_hot_in = prop_hot_in[u'x']
        if not u'q' in prop_hot_in: q_hot_in = 2
        else: q_hot_in = prop_hot_in[u'q']
        if not u'h' in prop_hot_in:
            multiRP.resetup(prop_hot_in, mRP=mRP)
        else: h_hot_in = prop_hot_in[u'h']

        #define max / min temp of fluid allowed for refprop
        #cold fluid (t_max)
        multiRP.resetup(prop_cold_in, mRP=mRP)
        if u'setmod' in prop_cold_in and u'htype' in prop_cold_in[u'setmod']:
            tmax_cold = (multiRP.limits(x_cold_in,
                prop_cold_in[u'setmod'][u'htype'])[u'tmax'] - dt)
        else: tmax_cold = multiRP.limits(x_cold_in)[u'tmax'] - dt

        #hot fluid (t_min)
        multiRP.resetup(prop_hot_in)
        if u'setmod' in prop_hot_in and u'htype' in prop_hot_in[u'setmod']:
            #check freezing point
            #IMPROVEMENT OPTION, refprop is incapable to calc. freezing points
            #at various mixtures, also the FLASH calculations are limited
            tmin_ht1 = max([multiRP.limitk(htype=prop_hot_in[u'setmod'][u'htype'],
                                            icomp=(i+1))[u'tmin'] \
                            for i in xrange(prop_hot_in[u'nc'])])
            #check limits
            tmin_ht2 = (multiRP.limits(x_hot_in,
                                       prop_hot_in[u'setmod'][u'htype'])[u'tmin'])
        else:
            #check freezing point
            #IMPROVEMENT OPTION, refprop is incapable to calc. freezing points
            #at various mixtures, also the FLASH calculations are limited
            tmin_ht1 = max([multiRP.limitk(icomp=(i+1))[u'tmin'] \
                            for i in xrange(prop_hot_in[u'nc'])])
            #check limits
            tmin_ht2 = multiRP.limits(x_hot_in)[u'tmin']
        tmin_hot = max(tmin_ht1, tmin_ht2) + dt

        #check if multiprocessing is allowed
        if SetMultiprocessing().__repr__() == u'on':
            ##start with multiprocessing##
            if mRP == None:
                #create if not exist
                mRP = multiRP.multirefprop()#initiate multirefprop

            #create children
            de1 = mRP[u'process'](target=_dE1, args=(Q_cold_in, p_cold_in,
                                                    dp_cold, t_hot_in, dt,
                                                    tmax_cold, x_cold_in,
                                                    h_cold_in, prop_cold_in,
                                                    mRP))
            de2 = mRP[u'process'](target=_dE2, args=(q_cold_in, Q_cold_in,
                                                    Q_hot_in, p_cold_in,
                                                    p_hot_in, x_cold_in,
                                                    x_hot_in, h_cold_in,
                                                    h_hot_in, tmin_hot,
                                                    prop_cold_in, prop_hot_in,
                                                    dp_cold, dp_hot, mRP))
            de4 = mRP[u'process'](target=_dE4, args=(dp_hot, dt, h_hot_in,
                                                    p_hot_in, prop_hot_in,
                                                    Q_hot_in, t_cold_in,
                                                    tmin_hot, x_hot_in, mRP))
            de6 = mRP[u'process'](target=_dE6, args=(h_cold_in, h_hot_in,
                                                    p_cold_in, p_hot_in,
                                                    prop_cold_in, prop_hot_in,
                                                    Q_cold_in,q_hot_in,
                                                    Q_hot_in, tmax_cold,
                                                    x_cold_in, x_hot_in,
                                                    dp_cold, dp_hot, mRP))
            de7 = mRP[u'process'](target=_dE7, args=(q_cold_in, Q_cold_in,
                                                    Q_hot_in, p_cold_in,
                                                    p_hot_in, x_cold_in,
                                                    x_hot_in, h_cold_in,
                                                    h_hot_in, tmin_hot,
                                                    prop_cold_in, prop_hot_in,
                                                    dp_cold, dp_hot, mRP))
            de8 = mRP[u'process'](target=_dE8, args=(h_cold_in, h_hot_in,
                                                    p_cold_in, p_hot_in,
                                                    prop_cold_in, prop_hot_in,
                                                    Q_cold_in, q_hot_in,
                                                    Q_hot_in, tmax_cold,
                                                    x_cold_in, x_hot_in,
                                                    dp_cold, dp_hot, mRP))
            dE_list = [de1, de2, de4, de6, de7, de8]

            #start and join children
            multiRP.run_mRP(dE_list)

            #assign results
            dE1 = mRP[u'result'][de1.name]
            dE2 = mRP[u'result'][de2.name]
            dE4 = mRP[u'result'][de4.name]
            dE6 = mRP[u'result'][de6.name]
            dE7 = mRP[u'result'][de7.name]
            dE8 = mRP[u'result'][de8.name]
            ##end with multiprocessing##

        #check if single processing is allowed
        elif SetMultiprocessing().__repr__() == u'off':
            dE1 = _dE1(Q_cold_in, p_cold_in, dp_cold, t_hot_in, dt, tmax_cold,
                       x_cold_in, h_cold_in, prop_cold_in, None)
            dE2 = _dE2(q_cold_in, Q_cold_in, Q_hot_in, p_cold_in, p_hot_in,
                       x_cold_in, x_hot_in, h_cold_in, h_hot_in, tmin_hot,
                       prop_cold_in, prop_hot_in, dp_cold, dp_hot, None)
            dE4 = _dE4(dp_hot, dt, h_hot_in, p_hot_in, prop_hot_in, Q_hot_in,
                       t_cold_in, tmin_hot, x_hot_in, None)
            dE6 = _dE6(h_cold_in, h_hot_in, p_cold_in, p_hot_in, prop_cold_in,
                       prop_hot_in, Q_cold_in,q_hot_in, Q_hot_in, tmax_cold,
                       x_cold_in, x_hot_in, dp_cold, dp_hot, None)
            dE7 = _dE7(q_cold_in, Q_cold_in, Q_hot_in, p_cold_in, p_hot_in,
                       x_cold_in, x_hot_in, h_cold_in, h_hot_in, tmin_hot,
                       prop_cold_in, prop_hot_in, dp_cold, dp_hot, None)
            dE8 = _dE8(h_cold_in, h_hot_in, p_cold_in, p_hot_in, prop_cold_in,
                       prop_hot_in, Q_cold_in, q_hot_in, Q_hot_in, tmax_cold,
                       x_cold_in, x_hot_in, dp_cold, dp_hot, None)


        #select the smallest dE
        dE = min(dE1, dE2, dE4, dE6, dE7, dE8)
        if dE < 0: dE = 0

        #calculate h_out
        #cold stream
        if prop_in[u't'] < prop_contra_in[u't']:
            h_out = h_in + (dE / Q_cold_in)
            multiRP.resetup(prop_cold_in)
            prop = multiRP.flsh(u'ph', p_out, h_out, x_cold_in)

        #hot stream
        elif prop_in[u't'] > prop_contra_in[u't']:
            h_out = h_in - (dE / Q_hot_in)
            multiRP.resetup(prop_hot_in)
            prop = multiRP.flsh(u'ph', p_out, h_out, x_hot_in)

    #restore Q value
    prop[u'Q'] = Q_in

    #set error margin
    if prop_in[u't'] < prop_contra_in[u't']:
        if prop[u'Q'] == scipy.inf:
            E_cold = scipy.inf
        else:
            E_cold = abs(prop_in[u'h'] - prop[u'h']) * prop[u'Q']
        if prop_contra_in[u'Q'] == scipy.inf:
            E_hot = scipy.inf
        else:
            E_hot = abs(prop_contra_in[u'h'] -
                        prop_contra_out[u'h']) * prop_contra_in[u'Q']
    else:
        if prop[u'Q'] == scipy.inf:
            E_hot = scipy.inf
        else:
            E_hot = abs(prop_in[u'h'] - prop[u'h']) * prop[u'Q']
        if prop_contra_in[u'Q'] == scipy.inf:
            E_cold = scipy.inf
        else:
            E_cold = abs(prop_contra_in[u'h'] -
                         prop_contra_out[u'h']) * prop_contra_in[u'Q']
    #correction to prevent div/0 error
    if E_cold == 0 or E_hot == 0:
        errvalue = 1
    elif E_cold == scipy.inf or E_hot == scipy.inf:
        errvalue = 0
    elif E_cold >= E_hot:
        errvalue = 1 - abs(E_hot / E_cold)
    else:
        errvalue = 1 - abs(E_cold / E_hot)

    equipspec ={u'type': u'EXCHANGER', u'name': name, u'dp': dp,
                u'dp_contra': dp_contra, u'dt': dt, u'fluidin': prop_in,
                u'fluidin_contra': prop_contra_in, u'fluidout': prop,
                u'fluidout_contra': prop_contra_out, u'errvalue': errvalue,
                u'E_hot':E_hot, u'E_cold':E_cold}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop

def _dE1(Q_cold_in, p_cold_in, dp_cold, t_hot_in, dt, tmax_cold, x_cold_in,
         h_cold_in, prop_cold_in, mRP):
    u"""define dE1 based on dt cold"""
    #resetup refprop
    multiRP.resetup(prop_cold_in, mRP=mRP)
    #check for infiniti
    if Q_cold_in == scipy.inf:
        dE1 = scipy.inf
        prop_dE1 = None
    #the calculations
    else:
        #determin pressure out
        p_cold_out = p_cold_in - dp_cold
        #pressure can not be negative
        if p_cold_out <= 0:
            p_cold_out = 0
        #determine temp. out
        t_cold_out = min(t_hot_in - dt, tmax_cold)
        #determine critical conditions
        crit = multiRP.critp(x_cold_in)
        if t_cold_out >= crit[u'tcrit'] or p_cold_out >= crit[u'pcrit']:
            #flash calculation
            h_cold_out = multiRP.flsh(u'tp', t_cold_out,
                                      p_cold_out,
                                      x_cold_in)[u'h']
        else:
            #liquid condition
            if t_cold_out <= multiRP.satp(p_cold_out, x_cold_in, 1)[u't']:
                D = multiRP.tprho(t_cold_out, p_cold_out, x_cold_in, 1)[u'D']
                h_cold_out = multiRP.enthal(t_cold_out, D, x_cold_in)[u'h']
            #vapor condition
            elif t_cold_out >= multiRP.satp(p_cold_out, x_cold_in, 2)[u't']:
                D = multiRP.tprho(t_cold_out, p_cold_out, x_cold_in, 2)[u'D']
                h_cold_out = multiRP.enthal(t_cold_out, D, x_cold_in)[u'h']
            #two phase condition
            else:
                try:
                    #flash calculation
                    h_cold_out = multiRP.flsh(u'tp', t_cold_out,
                                              p_cold_out, x_cold_in)[u'h']
                except:
                    #asume liquid condition
                    ############################################################
                    # -1e-05 added due to round off errors at sat. bounderies. #
                    ############################################################
                    if t_cold_out - 1e-05 <= multiRP.satp(p_cold_out,
                                                          x_cold_in, 1)[u't']:
                        D = multiRP.tprho(t_cold_out, p_cold_out,
                                          x_cold_in, 1)[u'D']
                        h_cold_out = multiRP.enthal(t_cold_out,
                                                    D, x_cold_in)[u'h']
                    else: raise EquipmentError()
        #heat exchanged
        dE1 = (h_cold_out - h_cold_in) * Q_cold_in
    #return results
    if SetMultiprocessing().__repr__() == u'on':
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = dE1
    elif SetMultiprocessing().__repr__() == u'off':
        return dE1

def _dE2(q_cold_in, Q_cold_in, Q_hot_in, p_cold_in, p_hot_in, x_cold_in,
         x_hot_in, h_cold_in, h_hot_in, tmin_hot, prop_cold_in, prop_hot_in,
         dp_cold, dp_hot, mRP):
    u"""define dE2 based on bubble point cold correction"""
    #settup cold fluid and initiate mRP
    multiRP.resetup(prop_cold_in, mRP=mRP)
    #check fluid not liquid
    phase = multiRP.getphase(prop_cold_in)
    if phase != u'liquid' and phase != u'saturated liquid':
        dE2 = scipy.inf
    #check infiniti
    elif Q_cold_in == scipy.inf or Q_hot_in == scipy.inf:
        dE2 = scipy.inf
    #calculate dE
    else:
        #determin pressure out
        p_cold_out = p_cold_in - dp_cold
        #pressure can not be negative
        if p_cold_out <= 0:
            p_cold_out = 0
        #p cold average
        p_cold_av = (p_cold_in + p_cold_out) / 2
        #calculate bubble point
        prop = multiRP.satp(p_cold_av, x_cold_in, 1)
        #temp at potential crossing of two fluids
        t_x_ing = prop[u't']
        #check for cross temperature below min temp.
        if t_x_ing < tmin_hot:
            dE2 = scipy.inf
        else:
            #calculate enthalpy at crossing (cold)
            h = multiRP.therm(t_x_ing, prop[u'Dliq'], x_cold_in)[u'h']
            #dE of cold part
            dE2a = (h - h_cold_in) * Q_cold_in
            #check for < 0 value
            if dE2a < 0: dE2 = scipy.inf
            else:
                #setup hot fluid
                multiRP.resetup(prop_hot_in, mRP=mRP)
                #determine critical conditions
                crit = multiRP.critp(x_hot_in)
                #determin pressure out
                p_hot_out = p_hot_in - dp_hot
                #pressure can not be negative
                if p_hot_out <= 0:
                    p_hot_out = 0
                #p cold average
                p_hot_av = (p_hot_in + p_hot_out) / 2
                if t_x_ing >= crit[u'tcrit'] or p_hot_av >= crit[u'pcrit']:
                    #flash calculation
                    h = multiRP.flsh(u'tp', t_x_ing,
                                     p_hot_av, x_hot_in)[u'h']
                else:
                    #liquid condition
                    if t_x_ing <= multiRP.satp(p_hot_av, x_hot_in, 1)[u't']:
                        D = multiRP.tprho(t_x_ing, p_hot_av, x_hot_in, 1)[u'D']
                        h = multiRP.enthal(t_x_ing, D, x_hot_in)[u'h']
                    #vapor condition
                    elif t_x_ing >= multiRP.satp(p_hot_av, x_hot_in, 2)[u't']:
                        D = multiRP.tprho(t_x_ing, p_hot_av, x_hot_in, 2)[u'D']
                        h = multiRP.enthal(t_x_ing, D, x_hot_in)[u'h']
                    #two phase condition
                    else:
                        #flash calculation
                        h = multiRP.flsh(u'tp', t_x_ing,
                                         p_hot_av, x_hot_in)[u'h']
                #dE of hot part
                dE2b = (h_hot_in - h) * Q_hot_in
                #check infiniti
                if dE2b < 0: dE2 = scipy.inf
                #sum dE hot and cold
                else: dE2 = dE2a + dE2b
    #return results
    if SetMultiprocessing().__repr__() == u'on':
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = dE2
    elif SetMultiprocessing().__repr__() == u'off':
        return dE2

def _dE4(dp_hot, dt, h_hot_in, p_hot_in, prop_hot_in, Q_hot_in, t_cold_in,
         tmin_hot, x_hot_in, mRP):
    u"""define dE4 based on dt hot"""
    #resetup and initiate mRP
    multiRP.resetup(prop_hot_in, mRP=mRP)
    #check for infiniti
    if Q_hot_in == scipy.inf:
        dE4 = scipy.inf
        prop_dE4 = None
    else:
        #calculate pressure out hot fluid
        p_hot_out = p_hot_in - dp_hot
        if p_hot_out <= 0:
            p_hot_out = 0
        #calculate temperature out hot fluid
        t_hot_out = max(t_cold_in + dt, tmin_hot)
        #determine critical conditions
        crit = multiRP.critp(x_hot_in)
        if t_hot_out >= crit[u'tcrit'] or p_hot_out >= crit[u'pcrit']:
            #flash calculation
            h_hot_out = multiRP.flsh(u'tp', t_hot_out, p_hot_out, x_hot_in)[u'h']
        else:
            #liquid condition
            ############################################################
            # -1e-07 added due to round off errors at sat. bounderies. #
            ############################################################
            if t_hot_out - 1e-07 <= multiRP.satp(p_hot_out, x_hot_in, 1)[u't']:
                D = multiRP.tprho(t_hot_out, p_hot_out, x_hot_in, 1)[u'D']
                h_hot_out = multiRP.enthal(t_hot_out, D, x_hot_in)[u'h']
            #vapor condition
            elif t_hot_out >= multiRP.satp(p_hot_out, x_hot_in, 2)[u't']:
                D = multiRP.tprho(t_hot_out, p_hot_out, x_hot_in, 2)[u'D']
                h_hot_out = multiRP.enthal(t_hot_out, D, x_hot_in)[u'h']
            #two phase condition
            else:
                try:
                    #flash calculation
                    h_hot_out = multiRP.flsh(u'tp', t_hot_out, p_hot_out,
                                             x_hot_in)[u'h']
                except:
                    #asume liquid condition
                    ############################################################
                    # -1e-05 added due to round off errors at sat. bounderies. #
                    ############################################################
                    if t_hot_out - 1e-05 <= multiRP.satp(p_hot_out, x_hot_in,
                                                         1)[u't']:
                        D = multiRP.tprho(t_hot_out, p_hot_out, x_hot_in,
                                          1)[u'D']
                        h_hot_out = multiRP.enthal(t_hot_out, D, x_hot_in)[u'h']
                    else: raise EquipmentError()
        #heat exchanged
        dE4 = (h_hot_in - h_hot_out) * Q_hot_in
    #return values
    if SetMultiprocessing().__repr__() == u'on':
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = dE4
    elif SetMultiprocessing().__repr__() == u'off':
        return dE4

def _dE6(h_cold_in, h_hot_in, p_cold_in, p_hot_in, prop_cold_in, prop_hot_in,
         Q_cold_in, q_hot_in, Q_hot_in, tmax_cold, x_cold_in, x_hot_in, dp_cold,
         dp_hot, mRP):
    u"""define dE6 based on dew point hot correction"""
    #setup hot fluid and initiate mRP
    multiRP.resetup(prop_hot_in, mRP=mRP)
    #check for q < 1
    phase = multiRP.getphase(prop_hot_in)
    if phase != u"vapor" and phase != u"saturated vapor" and phase != u"gas":
        dE6 = scipy.inf
    #check for infinity
    elif Q_cold_in == scipy.inf or Q_hot_in == scipy.inf:
        dE6 = scipy.inf
    else:
        #determin pressure out
        p_hot_out = p_hot_in - dp_hot
        #pressure can not be negative
        if p_hot_out <= 0:
            p_hot_out = 0
        #p cold average
        p_hot_av = (p_hot_in + p_hot_out) / 2
        #calculate bubble point
        prop = multiRP.satp(p_hot_av, x_hot_in, 2)
        #potential crossing temp
        t_x_ing = prop[u't']
        #check for temp > max allowed
        if t_x_ing > tmax_cold:
            dE6 = scipy.inf
        else:
            #enthalpy at crossing
            h = multiRP.therm(t_x_ing, prop[u'Dvap'], x_hot_in)[u'h']
            #dE hot fluid
            dE6a = (h_hot_in - h) * Q_hot_in
            #check for dE < 0
            if dE6a < 0: dE6 = scipy.inf
            else:
                #setup cold fluid
                multiRP.resetup(prop_cold_in)
                #determine critical conditions
                crit = multiRP.critp(x_cold_in)
                #determin pressure out
                p_cold_out = p_cold_in - dp_cold
                #pressure can not be negative
                if p_cold_out <= 0:
                    p_cold_out = 0
                #p cold average
                p_cold_av = (p_cold_in + p_cold_out) / 2
                if t_x_ing >= crit[u'tcrit'] or p_cold_av >= crit[u'pcrit']:
                    #flash calculation
                    h = multiRP.flsh(u'tp', t_x_ing, p_cold_av, x_cold_in)[u'h']
                else:
                    #liquid condition
                    if t_x_ing <= multiRP.satp(p_cold_av, x_cold_in, 1)[u't']:
                        D = multiRP.tprho(t_x_ing, p_cold_av, x_cold_in, 1)[u'D']
                        h = multiRP.enthal(t_x_ing, D, x_cold_in)[u'h']
                    #vapor condition
                    elif t_x_ing >= multiRP.satp(p_cold_av, x_cold_in, 2)[u't']:
                        D = multiRP.tprho(t_x_ing, p_cold_av, x_cold_in, 2)[u'D']
                        h = multiRP.enthal(t_x_ing, D, x_cold_in)[u'h']
                    #two phase condition
                    else:
                        #flash calculation
                        h = multiRP.flsh(u'tp', t_x_ing, p_cold_av,
                                         x_cold_in)[u'h']
                #dE cold fluid
                dE6b = (h - h_cold_in) * Q_cold_in
                #check for dE < 0
                if dE6b < 0: dE6 = scipy.inf
                #calculate dE
                else: dE6 = dE6a + dE6b
    #return results
    if SetMultiprocessing().__repr__() == u'on':
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = dE6
    elif SetMultiprocessing().__repr__() == u'off':
        return dE6

def _dE7(q_cold_in, Q_cold_in, Q_hot_in, p_cold_in, p_hot_in, x_cold_in,
         x_hot_in, h_cold_in, h_hot_in, tmin_hot, prop_cold_in, prop_hot_in,
         dp_cold, dp_hot, mRP):
    u"""define dE7 based on critical point cold correction"""
    #settup cold fluid and initiate mRP
    multiRP.resetup(prop_cold_in, mRP=mRP)
    #check fluid above pcrit.
    phase = multiRP.getphase(prop_cold_in)
    if phase != u"compressible liquid":
        dE7 = scipy.inf
    #check infiniti
    elif Q_cold_in == scipy.inf or Q_hot_in == scipy.inf:
        dE7 = scipy.inf
    #calculate dE
    else:
        #calculate critical cond.
        prop = multiRP.critp(x_cold_in)
        #temp at potential crossing of two fluids
        t_x_ing = prop[u'tcrit']
        #check for cross temperature below min temp.
        if t_x_ing < tmin_hot:
            dE7 = scipy.inf
        else:
            #determin pressure out
            p_cold_out = p_cold_in - dp_cold
            #pressure can not be negative
            if p_cold_out <= 0:
                p_cold_out = 0
            #p cold average
            p_cold_av = (p_cold_in + p_cold_out) / 2
            #calculate enthalpy at crossing (cold)
            h = multiRP.flsh(u'tp', t_x_ing, p_cold_av, x_cold_in)[u'h']
            #dE of cold part
            dE7a = (h - h_cold_in) * Q_cold_in
            #check for < 0 value
            if dE7a < 0: dE7 = scipy.inf
            else:
                #setup hot fluid
                multiRP.resetup(prop_hot_in, mRP=mRP)
                #critical cond.
                crit = multiRP.critp(x_hot_in)
                #determin pressure out
                p_hot_out = p_hot_in - dp_hot
                #pressure can not be negative
                if p_hot_out <= 0:
                    p_hot_out = 0
                #p cold average
                p_hot_av = (p_hot_in + p_hot_out) / 2
                if t_x_ing >= crit[u'tcrit'] or p_hot_av >= crit[u'pcrit']:
                    #flash calculation
                    h = multiRP.flsh(u'tp', t_x_ing, p_hot_av, x_hot_in)[u'h']
                else:
                    #liquid condition
                    if t_x_ing <= multiRP.satp(p_hot_av, x_hot_in, 1)[u't']:
                        D = multiRP.tprho(t_x_ing, p_hot_av, x_hot_in, 1)[u'D']
                        h = multiRP.enthal(t_x_ing, D, x_hot_in)[u'h']
                    #vapor condition
                    elif t_x_ing >= multiRP.satp(p_hot_av, x_hot_in, 2)[u't']:
                        D = multiRP.tprho(t_x_ing, p_hot_av, x_hot_in, 2)[u'D']
                        h = multiRP.enthal(t_x_ing, D, x_hot_in)[u'h']
                    #two phase condition
                    else:
                        #flash calculation
                        h = multiRP.flsh(u'tp', t_x_ing, p_hot_av, x_hot_in)[u'h']
                #dE of hot part
                dE7b = (h_hot_in - h) * Q_hot_in
                #check infiniti
                if dE7b < 0: dE7 = scipy.inf
                #sum dE hot and cold
                else: dE7 = dE7a + dE7b
    #return results
    if SetMultiprocessing().__repr__() == u'on':
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = dE7
    elif SetMultiprocessing().__repr__() == u'off':
        return dE7

def _dE8(h_cold_in, h_hot_in, p_cold_in, p_hot_in, prop_cold_in, prop_hot_in,
         Q_cold_in, q_hot_in, Q_hot_in, tmax_cold, x_cold_in, x_hot_in, dp_cold,
         dp_hot, mRP):
    u"""define dE8 based on critical point hot correction"""
    #setup hot fluid and initiate mRP
    multiRP.resetup(prop_hot_in, mRP=mRP)
    #check for q < 1
    phase = multiRP.getphase(prop_hot_in)
    if phase != u"Supercritical fluid":
        dE8 = scipy.inf
    #check for infinity
    elif Q_cold_in == scipy.inf or Q_hot_in == scipy.inf:
        dE8 = scipy.inf
    else:
        #calculate critical point
        prop = multiRP.critp(x_hot_in)
        #potential crossing temp
        t_x_ing = prop[u'tcrit']
        #check for temp > max allowed
        if t_x_ing > tmax_cold:
            dE6 = scipy.inf
        else:
            #determin pressure out
            p_hot_out = p_hot_in - dp_hot
            #pressure can not be negative
            if p_hot_out <= 0:
                p_hot_out = 0
            #p cold average
            p_hot_av = (p_hot_in + p_hot_out) / 2
            #enthalpy at crossing
            h = multiRPflsh(u'tp', t_x_ing, p_hot_av, x_hot_in)[u'h']
            #dE hot fluid
            dE8a = (h_hot_in - h) * Q_hot_in
            #check for dE < 0
            if dE8a < 0: dE8 = scipy.inf
            else:
                #setup cold fluid
                multiRP.resetup(prop_cold_in)
                #determine critical conditions
                crit = multiRP.critp(x_cold_in)
                #determin pressure out
                p_cold_out = p_cold_in - dp_cold
                #pressure can not be negative
                if p_cold_out <= 0:
                    p_cold_out = 0
                #p cold average
                p_cold_av = (p_cold_in + p_cold_out) / 2
                if t_x_ing >= crit[u'tcrit'] or p_cold_av >= crit[u'pcrit']:
                    #flash calculation
                    h = multiRP.flsh(u'tp', t_x_ing, p_cold_av, x_cold_in)[u'h']
                else:
                    #liquid condition
                    if t_x_ing <= multiRP.satp(p_cold_av, x_cold_in, 1)[u't']:
                        D = multiRP.tprho(t_x_ing, p_cold_av, x_cold_in, 1)[u'D']
                        h = multiRP.enthal(t_x_ing, D, x_cold_in)[u'h']
                    #vapor condition
                    elif t_x_ing >= multiRP.satp(p_cold_av, x_cold_in, 2)[u't']:
                        D = multiRP.tprho(t_x_ing, p_cold_av, x_cold_in, 2)[u'D']
                        h = multiRP.enthal(t_x_ing, D, x_cold_in)[u'h']
                    #two phase condition
                    else:
                        #flash calculation
                        h = multiRP.flsh(u'tp', t_x_ing, p_cold_av,
                                         x_cold_in)[u'h']
                #dE cold fluid
                dE8b = (h - h_cold_in) * Q_cold_in
                #check for dE < 0
                if dE8b < 0: dE8 = scipy.inf
                #calculate dE
                else: dE8 = dE8a + dE8b
    #return results
    if SetMultiprocessing().__repr__() == u'on':
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = dE8
    elif SetMultiprocessing().__repr__() == u'off':
        return dE8


def turbine_separator(prop_in, p_out, dp, eff=1.0, e_eff=1.0, wet=0.0,
                      name=False, mRP=None):
    u"""Calculate fluid properties at turbine_separator outlet
    (liquid and vapor).
    equipment specification out put through TurbineSeparator

    input:
        prop_in--fluid properties library contain keys:
            'Q'--mole flow (mole / sec)
            'p'--pressure
            'x'--mole fraction
            'h'--enthalpy
        p_out--turbine output pressure
        dp--pressure differential across separator
        eff--efficiency of turbine (std = 1)
        e_eff--efficiency of alternator (std = 1)
        wet--Turbine max. water content in outlet
    output:
        prop--fluid properties
        equipspec--equipment details"""
    _inputerrorcheck(locals(), u'turbine_separator')
    multiRP.resetup(prop_in, mRP=mRP)

    #separator
    prop_liq, prop_vap = separator(prop_in, dp, mRP=mRP)

    ##start with multiprocessing##
    if SetMultiprocessing().__repr__() == u'on':
        if mRP == None:
            #create if not exist
            mRP = multiRP.multirefprop()#initiate multirefprop

        #create children
        tursep_reg = mRP[u'process'](target=reg_valve, args=(prop_liq, p_out),
                                    kwargs={u'mRP':mRP})
        tursep_tur = mRP[u'process'](target=turbine, args=(prop_vap, p_out),
                                    kwargs={u'mRP':mRP, u'eff':eff, u'e_eff':e_eff,
                                            u'wet':wet})
        tursep_list = [tursep_reg, tursep_tur]

        #start and join children
        multiRP.run_mRP(tursep_list)

        #assign results
        prop_reg = mRP[u'result'][tursep_reg.name]
        prop_tur = mRP[u'result'][tursep_tur.name]
        ##end with multiprocessing##

    #or if single core
    elif SetMultiprocessing().__repr__() == u'off':
        prop_reg = reg_valve(prop_liq, p_out, mRP=None)
        prop_tur = turbine(prop_vap, p_out, mRP=None, eff=eff, e_eff=e_eff,
                           wet=wet)

    #set increased back pressure
    if u'backpress' in prop_tur:
        incr_back_press = prop_tur[u'backpress']
        prop_tur.pop(u'backpress')

    #flowmerge
    if prop_tur[u'Q'] == 0:
        prop = prop_reg
    elif prop_reg[u'Q'] == 0:
        prop = prop_tur
    else:
        prop = flowmerge(prop_reg, prop_tur)

    #power generated
    pwr = ((prop_vap[u'h'] - prop_tur[u'h']) * prop_vap[u'Q']) * e_eff

    #create equipment spec
    equipspec ={u'type': u'TURBINESEP', u'name': name, u'eff': eff, u'p_out': p_out,
                u'fluidin': prop_in, u'fluidout': prop, u'pwr': pwr, u'dp': dp,
                u'e_eff':e_eff, u'wet':wet, u'incr_back_press':incr_back_press}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop


def turbine_bleed_sep(prop_in, p_out, P_outbleed, dp, Q_ratio, Q_ratiobleed,
                      eff=1.0, e_eff=1.0, wet=0.0, name=False, mRP=None):
    u"""Calculate fluid properties at turbine_separator outlet
    (liquid and vapor).
    equipment specification out put through TurbineSeparator

    input:
        prop_in--fluid properties library contain keys:
            'Q'--mole flow (mole / sec)
            'p'--pressure
            'x'--mole fraction
            'h'--enthalpy
        p_out--turbine output pressure
        p_outbleed--turbine output pressure at bleed
        dp--pressure differential across separator
        Q_ratio--output flow of bleed as a ratio of the total gasflow
            0 = 0% bleed flow (100% main flow)
            1 = 100% bleed flow (0% main flow)
        Q_ratiobleed--condensate flow ratio (from separator)
            0 = 0% condensate to bleed flow (100% to main flow)
            1 = 100% condensate to bleed flow (0% to main flow)
        eff--efficiency of turbine (std = 1)
        e_eff--efficiency of alternator (std = 1)
        wet--Turbine max. water content in outlet
    output:
        prop--fluid properties main turbine
        prop_bleed--fluid properties at bleed out (incl. sep. liquid)
        equipspec--equipment details"""
    _inputerrorcheck(locals(), u'turbine_bleed_sep')
    multiRP.resetup(prop_in, mRP=mRP)

    #separator
    prop_liq, prop_vap = separator(prop_in, dp)

    #check if fluid in only liquid
    if prop_vap[u'Q'] == 0:
        ##start with multiprocessing##
        if SetMultiprocessing().__repr__() == u'on':
            if mRP == None:
                #create if not exist
                mRP = multiRP.multirefprop()#initiate multirefprop

            #create children
            tursep_0vap = mRP[u'process'](target=reg_valve,
                                         args=(prop_liq, p_out),
                                         kwargs={u'mRP':mRP})
            tursep_0vapbleed = mRP[u'process'](target=reg_valve,
                                              args=(prop_liq, p_outbleed),
                                              kwargs={u'mRP':mRP})
            tursep_list = [tursep_0vap, tursep_0vapbleed]

            #start and join children
            multiRP.run_mRP(tursep_list)

            #assign results
            prop = mRP[u'result'][tursep_0vap.name]
            prop_bleed = mRP[u'result'][tursep_0vapbleed.name]
            ##end with multiprocessing##

        #single core
        elif SetMultiprocessing().__repr__() == u'off':
            prop = reg_valve(prop_liq, p_out, mRP=None)
            prop_bleed = reg_valve(prop_liq, p_outbleed, mRP=None)

        #correct flow value
        prop[u'Q'] *= 1 - Q_ratiobleed
        prop_bleed[u'Q'] *= Q_ratiobleed

    #check if fluid contains vapor
    else:
        #copy prop_vap to bleed properties and modify Q value
        prop_bleed = copy.copy(prop_vap)
        prop_vap[u'Q'] *= 1 - Q_ratio
        prop_bleed[u'Q'] *= Q_ratio
        #copy prop_liq to bleed properties and modify Q value
        prop_liq_bleed = copy.copy(prop_liq)
        prop_liq[u'Q'] *= 1 - Q_ratiobleed
        prop_liq_bleed[u'Q'] *= Q_ratiobleed

        ##start with multiprocessing##
        if SetMultiprocessing().__repr__() == u'on':
            if mRP == None:
                #create if not exist
                mRP = multiRP.multirefprop()#initiate multirefprop

            #create children
            tursep_reg = mRP[u'process'](target=reg_valve,
                                        args=(prop_liq, p_out),
                                        kwargs={u'mRP':mRP})
            tursep_reg_bleed = mRP[u'process'](target=reg_valve,
                                              args=(prop_liq_bleed, p_outbleed),
                                              kwargs={u'mRP':mRP})
            tursep_tur = mRP[u'process'](target=turbine,
                                        args=(prop_vap, p_out),
                                        kwargs={u'mRP':mRP, u'eff':eff,
                                                u'e_eff':e_eff, u'wet':wet})
            tursep_bleed = mRP[u'process'](target=turbine,
                                          args=(prop_bleed, p_outbleed),
                                          kwargs={u'mRP':mRP, u'eff':eff,
                                                  u'e_eff':e_eff, u'wet':wet})
            tursep_list = [tursep_reg, tursep_reg_bleed, tursep_tur,
                           tursep_bleed]

            #start and join children
            multiRP.run_mRP(tursep_list)

            #assign results
            prop_reg = mRP[u'result'][tursep_reg.name]
            prop_reg_bleed = mRP[u'result'][tursep_reg_bleed.name]
            prop_tur = mRP[u'result'][tursep_tur.name]
            prop_bleed = mRP[u'result'][tursep_bleed.name]
            ##end with multiprocessing##

        #single core
        elif SetMultiprocessing().__repr__() == u'off':
            prop_reg = reg_valve(prop_liq, p_out, mRP=None)
            prop_reg_bleed = reg_valve(prop_liq_bleed, p_outbleed, mRP=None)
            prop_tur = turbine(prop_vap, p_out, mRP=None, eff=eff, e_eff=e_eff,
                               wet=wet)
            prop_bleed = turbine(prop_bleed, p_outbleed, mRP=None, eff=eff,
                                 e_eff=e_eff, wet=wet)

        #set increased back pressure
        incr_back_press = prop_tur[u'backpress']
        incr_back_press_bleed = prop_bleed[u'backpress']

        #flowmerge
        ##start with multiprocessing##
        if SetMultiprocessing().__repr__() == u'on':
            if mRP == None:
                #create if not exist
                mRP = multiRP.multirefprop()#initiate multirefprop

            #create children
            tursep_tur = mRP[u'process'](target=flowmerge,
                                        args=(prop_reg, prop_tur),
                                        kwargs={u'mRP':mRP})
            tursep_bleed = mRP[u'process'](target=flowmerge,
                                          args=(prop_reg_bleed, prop_bleed),
                                          kwargs={u'mRP':mRP})
            tursep_list = [tursep_tur, tursep_bleed]

            #start and join children
            multiRP.run_mRP(tursep_list)

            #assign results
            prop = mRP[u'result'][tursep_tur.name]
            prop_bleed = mRP[u'result'][tursep_bleed.name]
            ##end multiprocessing##

        #single core
        elif SetMultiprocessing().__repr__() == u'off':
            prop = flowmerge(prop_reg, prop_tur, mRP=None)
            prop_bleed = flowmerge(prop_reg_bleed, prop_bleed, mRP=None)

    #power generated
    pwr = (((prop_vap[u'h'] - prop_tur[u'h']) * prop_vap[u'Q']) * e_eff +
          ((prop_vap[u'h'] - prop_bleed[u'h']) * prop_bleed[u'Q']) * e_eff)

    #create equipment spec
    equipspec ={u'type': u'TURBINEBLEEDSEP', u'name': name, u'eff': eff,
                u'p_out': p_out, u'fluidin': prop_in, u'fluidout': prop,
                u'pwr': pwr, u'dp': dp, u'e_eff':e_eff, u'p_outbleed':p_outbleed,
                u'Q_ratio':Q_ratio, u'fluidbleedout':prop_bleed,
                u'Q_ratiobleed':Q_ratiobleed, u'wet':wet,
                u'incr_back_press':incr_back_press,
                u'incr_back_press_bleed':incr_back_press_bleed}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = (prop,
                                                                    prop_bleed)

    #return values
    if name != False:
        return prop, prop_bleed, equipspec
    else: return prop, prop_bleed


def reg_valve(prop_in, p_out, name=False, mRP=None):
    u"""Calculate fluid properties downstream press. reg. valve with constant
    enthalpy (100% efficency).
    Equipment specification output through RegValve

    input:
        prop_in--fluid properties library contain keys:
            'h'--enthalpy
            'x'--fluid composition
        p_out--set pressure downstream the press. reg. valve
    output:
        prop--fluid properties
        equipspec--equipment details"""
    _inputerrorcheck(locals(), u'reg_valve')
    multiRP.resetup(prop_in, mRP=mRP)
    h = prop_in[u'h']
    x = prop_in[u'x']

    #calculate properties downstream valve
    prop = multiRP.flsh(u'ph', p_out, h, x)

    #add Q
    if u'Q' in prop_in:
        prop[u'Q'] = prop_in[u'Q']

    equipspec ={u'type': u'REGVALVE', u'name': name, u'p_out': p_out,
                u'fluidin': prop_in, u'fluidout': prop}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop


def separator(prop_in, dp, name=False, mRP=None):
    u"""Calculate fluid properties at separator outlets (liquid and vapor).
    equipment specification out put through Separator

    input:
        prop_in--fluid properties library contain keys:
            'Q'--mole flow (mole / sec)
            'p'--pressure
            'x'--mole fraction
            'h'--enthalpy
        dp--pressure differential across separator
    output:
        prop_liq--Liquid fluid properties
        prop_vap--vapor fluid properties
        equipspec--equipment details"""
    _inputerrorcheck(locals(), u'separator')
    multiRP.resetup(prop_in, mRP=mRP)
    p = prop_in[u'p'] - dp
    x = prop_in[u'x']
    h = prop_in[u'h']
    Q = prop_in[u'Q']

    #calculate properties with corrected p and constant h
    #enhance speed
    if dp > 0:
        #inhibit error reporting
        if unicode(multiRP.SetErrorDebug()) == u'on':
            sed = multiRP.SetErrorDebug.on
            multiRP.SetErrorDebug.off()
        elif unicode(multiRP.SetErrorDebug()) == u'off':
            sed = multiRP.SetErrorDebug.off
        try:
            prop = multiRP.ph2ph(p, h, x)
            prop.update(sed())
        except multiRP.RefpropError:
            sed() #re-activate error reporting
            prop = multiRP.flsh(u'ph', p, h, x)
    else: prop = copy.copy(prop_in)

    prop_liq = copy.copy(prop)
    prop_vap = copy.copy(prop)
    prop_liq.update({u'D': prop[u'Dliq'], u'x': prop[u'xliq']})
    prop_vap.update({u'D': prop[u'Dvap'], u'x': prop[u'xvap']})

    #vapor quality to be between 0 and 1
    phase = multiRP.getphase(prop)
    if phase == u"liquid" or phase == u"compressible liquid" \
    or phase == u"saturated liquid":
        q = 0
    elif phase == u"saturated vapor" or phase == u"gas" \
    or phase == u"Supercritical fluid" or phase == u"vapor":
        q = 1
    elif phase == u"2 phase":
        q = prop[u'q']

    ##start with multiprocessing##
    if SetMultiprocessing().__repr__() == u'on':
        if mRP == None:
            #create mRP if not exist
            mRP = multiRP.multirefprop()#initiate multirefprop

        #create children
        sep_liq = mRP[u'process'](target=_separator_liq, args=(prop_liq, p, mRP))
        sep_vap = mRP[u'process'](target=_separator_vap, args=(prop_vap, p, mRP))
        sep_list = [sep_liq, sep_vap]

        #start and join children
        multiRP.run_mRP(sep_list)

        #assign results
        prop_liq = mRP[u'result'][sep_liq.name]
        prop_vap = mRP[u'result'][sep_vap.name]
        ##end with multiprocessing##

    #single core
    elif SetMultiprocessing().__repr__() == u'off':
        prop_liq = _separator_liq(prop_liq, p, None)
        prop_vap = _separator_vap(prop_vap, p, None)

    #add Q to input values for liquid and vapor
    prop_liq[u'Q'] = Q * (1 - q)
    prop_vap[u'Q'] = Q * q

    #create equipspec
    equipspec ={u'type': u'SEPARATOR', u'name': name, u'dp': dp,
                u'fluidin': prop_in, u'fluidliq': prop_liq, u'fluidvap': prop_vap}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = prop

    #return values
    if name != False:
        return prop_liq, prop_vap, equipspec
    else: return prop_liq, prop_vap

#child1 for separator, calculate liquid properties
def _separator_liq(prop_liq, p, mRP):
    u'child calculation for def "separator"'
    prop = multiRP.therm(prop_liq[u't'], prop_liq[u'D'], prop_liq[u'x'],
                         prop=prop_liq, mRP=mRP)
    prop[u'p'] = p
    prop[u'q'] = 0 #to saturated liquid state
    if SetMultiprocessing().__repr__() == u'on':
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = prop
    elif SetMultiprocessing().__repr__() == u'off':
        return prop

#child2 for separator, calculate vapor properties
def _separator_vap(prop_vap, p, mRP):
    u'child calculation for def "separator"'
    prop = multiRP.therm(prop_vap[u't'], prop_vap[u'D'], prop_vap[u'x'],
                         prop=prop_vap, mRP=mRP)
    prop[u'p'] = p
    prop[u'q'] = +1 #to saturated vapor state
    if SetMultiprocessing().__repr__() == u'on':
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = prop
    elif SetMultiprocessing().__repr__() == u'off':
        return prop


def compressor(prop_in, dp, eff=1.0, e_eff=1.0, name=False, mRP=None):
    u"""Calculate fluid properties at compressor discharge.
    Equipment specification output through Compressor

    input:
        prop_in--fluid properties library contain keys:
            'Q'--mole flow (mole / sec)
            's'--entropy
            'p'--pressure
            'x'--mole fraction
            'h'--enthalpy (optional)
        dp--compressor head
        eff--pump efficiency (std = 1)
        e_eff--motor efficiency (std = 1)
    output:
        prop--fluid properties
        equipspec--equipment details"""
    _inputerrorcheck(locals(), u'compressor')
    multiRP.resetup(prop_in, mRP=mRP)
    p = prop_in[u'p'] + dp
    s = prop_in[u's']
    x = prop_in[u'x']
    flshcalc = False
    if u'h' in prop_in: h = prop_in[u'h']
    else:
        #inhibit error reporting
        if unicode(multiRP.SetErrorDebug()) == u'on':
            sed = multiRP.SetErrorDebug.on
            multiRP.SetErrorDebug.off()
        elif unicode(multiRP.SetErrorDebug()) == u'off':
            sed = multiRP.SetErrorDebug.off
        try: #use flsh1 calculation for speed increase
            h = prop = multiRP.psvap(prop_in[u'p'], s, x)[u'h']
            prop.update(sed())
        except multiRP.RefpropError:
            sed() #re-activate error reporting
            flshcalc = True
            h = multiRP.flsh(u'ps', p, s, x)[u'h']
    Q = prop_in[u'Q']

    #inhibit error reporting
    if unicode(multiRP.SetErrorDebug()) == u'on':
        sed = multiRP.SetErrorDebug.on
        multiRP.SetErrorDebug.off()
    elif unicode(multiRP.SetErrorDebug()) == u'off':
        sed = multiRP.SetErrorDebug.off
    #calculate enthalpy at 100% efficiency
    try: #use flsh1 calculation for speed increase
        if flshcalc: raise multiRP.RefpropError
        prop = multiRP.psvap(p, s, x)
        prop.update(sed())
    except multiRP.RefpropError:
        sed() #re-activate error reporting
        flshcalc = True
        prop = multiRP.flsh(u'ps', p, s, x)

    #correct h for efficiency losses
    h_out = h + ((prop[u'h'] - h) / eff)

    #calculate properties with corrected h and discharge p
    #inhibit error reporting
    if unicode(multiRP.SetErrorDebug()) == u'on':
        sed = multiRP.SetErrorDebug.on
        multiRP.SetErrorDebug.off()
    elif unicode(multiRP.SetErrorDebug()) == u'off':
        sed = multiRP.SetErrorDebug.off
    try: #use flsh1 calculation for speed increase
        if flshcalc: raise multiRP.RefpropError
        prop = multiRP.phvap(p, h_out, x)
        prop.update(sed())
    except multiRP.RefpropError:
        sed() #re-activate error reporting
        prop = multiRP.flsh(u'ph', p, h_out, x)

    #restore Q value
    prop[u'Q'] = Q

    #power consumption
    pwr = ((h_out - h) * Q) / e_eff

    equipspec ={u'type': u'COMPRESSOR', u'name': name, u'dp': dp, u'eff': eff,
                u'fluidin': prop_in, u'fluidout': prop, u'pwr': -pwr,
                u'e_eff':e_eff}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop



def subcool(prop_in, dt, name=False, mRP=None):
    u"""Calculate fluid properties based on subcooled temperature

    input:
        prop--fluid properties library contain keys:
            't'--subcool temperature.
            'x'--fluid composition
        dt--temperature difference between bubble point and subcool
    output:
        prop--fluid out properties
        equipspec--equipment details"""
    _inputerrorcheck(locals(), u'subcool')
    multiRP.resetup(prop_in, mRP=mRP)

    #calculate saturated properties at elevated temperature
    p = multiRP.satt(prop_in[u't'] + dt, prop_in[u'x'], 1)[u'p']

    #calculate sub cooled properties at set temp and calc. press
    try:
        D = multiRP.tprho(prop_in[u't'], p, prop_in[u'x'], 1)[u'D']
        prop = multiRP.therm(prop_in[u't'], D, prop_in[u'x'])
        prop[u'q'] = -1
        prop[u'p'] = p
    except:
        prop = multiRP.flsh(u'tp', prop_in[u't'], p, prop_in[u'x'])

    #add Q
    if u'Q' in prop_in:
        prop[u'Q'] = prop_in[u'Q']

    #create equipment spec
    equipspec ={u'type': u'SUBCOOL', u'name': name, u'dt': dt, u'fluidin': prop_in,
                u'fluidout': prop}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop


def turbine(prop_in, p_out, eff=1.0, e_eff=1.0, wet=0.0, name=False, mRP=None):
    u"""Calculate fluid properties at turbine discharge.
    equipment specification out put through Turbine

    input:
        prop_in--fluid properties library contain keys:
            'Q'--mole flow (mole / sec)
            's'--entropy
            'x'--mole fraction
            'h'--enthalpy
        p_out--turbine outlet pressure
        eff--pump efficiency (std = 1)
        e_eff--Alternator efficiency (std = 1)
        wet--Turbine max. water content in outlet
    output:
        prop_in--fluid properties
        equipspec--equipment details"""
    _inputerrorcheck(locals(), u'turbine')
    multiRP.resetup(prop_in, mRP=mRP)
    p = p_out
    s = prop_in[u's']
    x = prop_in[u'x']
    h = prop_in[u'h']
    Q = prop_in[u'Q']

    def _turbine():
        #set error
        flsh_2ph = False
        flsh_vap = False
        #calculate enthalpy at 100% efficiency
        #inhibit error reporting
        if unicode(multiRP.SetErrorDebug()) == u'on':
            sed = multiRP.SetErrorDebug.on
            multiRP.SetErrorDebug.off()
        elif unicode(multiRP.SetErrorDebug()) == u'off':
            sed = multiRP.SetErrorDebug.off
        #try flash1 calculation
        try:
            prop = multiRP.ps2ph(p, s, x)
            prop.update(sed())
        except multiRP.RefpropError:
            flsh_2ph = True
            try:
                prop = multiRP.psvap(p, s, x)
                prop.update(sed())
            except multiRP.RefpropError:
                sed()
                flsh_vap = True
                prop = multiRP.flsh(u'ps', p, s, x)
        #correct h for efficiency losses
        h_out = h - ((h - prop[u'h']) * eff)

        #calculate properties with corrected h and discharge p
        #inhibit error reporting
        if unicode(multiRP.SetErrorDebug()) == u'on':
            sed = multiRP.SetErrorDebug.on
            multiRP.SetErrorDebug.off()
        elif unicode(multiRP.SetErrorDebug()) == u'off':
            sed = multiRP.SetErrorDebug.off
        try:
            if flsh_2ph: raise multiRP.RefpropError
            prop = multiRP.ph2ph(p, h_out, x)
            prop.update(sed())
        except multiRP.RefpropError:
            try:
                if flsh_vap: raise multiRP.RefpropError
                prop = multiRP.phvap(p, h_out, x)
                prop.update(sed())
            except multiRP.RefpropError:
                sed()
                prop = multiRP.flsh(u'ph', p, h_out, x)
        #restore input values
        prop[u'Q'] = Q
        #power generated
        pwr = ((h - h_out) * prop_in[u'Q']) * e_eff
        #return
        return (prop, pwr)

    #run turbine
    prop, pwr = _turbine()

    #check wet gas within limits
    if prop[u'q'] >= wet:
        incr_back_press = 0
    else:
        #set initial press adjustment value
        p_adj = 1
        #set initial delta pressure
        press_delta = prop_in[u'p'] - p
        #create while loop to find optimum back pressure
        while not wet - 0.00001 < prop[u'q'] < wet + 0.00001:
            #adjust press_adjust with factor 2
            p_adj *= 2
            #store previous run
            q_old = prop[u'q']
            p_old = p
            #correct pressure delta value
            if prop[u'q'] <= wet:
                press_delta -= press_delta / p_adj
            elif prop[u'q'] >= wet:
                press_delta += press_delta / p_adj
            #new pressure value
            p = prop_in[u'p'] - press_delta
            #run turbine
            prop, pwr = _turbine()
            #compare new values against old to determin insoluble situation
            if q_old > prop[u'q'] and p_old < p:
                raise EquipmentinfeasibleError(u'error raised insoluble q value')
            if p_adj >= 2 ** 16:
                raise EquipmentinfeasibleError(u'error raised unsolvable q value')
        #additional press
        incr_back_press = p - p_out
        #depress to reach p_out
        prop = reg_valve(prop, p_out, mRP=mRP)

    #create equipment spec
    equipspec ={u'type': u'TURBINE', u'name': name, u'eff': eff, u'p_out': p_out,
                u'fluidin': prop_in, u'fluidout': prop, u'pwr': pwr,
                u'e_eff':e_eff, u'wet':wet, u'incr_back_press':incr_back_press}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        prop[u'backpress'] = incr_back_press
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = prop
        prop.pop(u'backpress')

    #return values
    if name != False:
        return prop, equipspec
    else:
        #add back pressure for turbine separator
        prop[u'backpress'] = incr_back_press
        return prop


def flowmerge(prop_1, prop_2=None, name=False, mRP=None):
    u"""Calculate fluid properties of 2 flow merge.
    equipment specification out put through FlowMerge

    input:
        prop_1--fluid properties library contain keys:
            'Q'--mole flow (mole / sec)
            'p'--pressure
            'x'--mole fraction
            'h'--enthalpy (optional)
        prop_2--fluid properties library contain keys:
            'Q'--mole flow (mole / sec)
            'p'--pressure
            'x'--mole fraction
            'h'--enthalpy (optional)
    output:
        prop--fluid properties
        equipspec--equipment details"""
    #del prop_2 if None to avoid inputerrorcheck
    if prop_2 == None:
        del prop_2
    _inputerrorcheck(locals(), u'flowmerge')

    #return prop_1 value if prop_2 had been deleted
    if not u'prop_2' in locals():
        prop = copy.copy(prop_1)

        #create equipment spec
        equipspec ={u'type': u'FLOWMERGE', u'name': name, u'fluid1': prop_1,
                    u'fluid2': None, u'fluidout': prop}

    else:
        #confirm both prop input are same refprop setup
        if multiRP.setup_details(prop_1) != multiRP.setup_details(prop_2):
            raise EquipmentinputError(u'''prop_1 and prop_2 input are not
                                      matching, ensure input are of same
                                      setup''')

        multiRP.resetup(prop_1, mRP=mRP)
        p1 = prop_1[u'p']
        x1 = prop_1[u'x']
        h1 = decimal.Decimal(prop_1[u'h'])
        Q1 = decimal.Decimal(prop_1[u'Q'])
        p2 = prop_2[u'p']
        x2 = prop_2[u'x']
        h2 = decimal.Decimal(prop_2[u'h'])
        Q2 = decimal.Decimal(prop_2[u'Q'])
        Q = Q1 + Q2

        #calculate merged x value.
        x = [(x1[each] * Q1 + x2[each] * Q2) / Q for each in xrange(len(x1))]
        x = multiRP.normalize(x)[u'x']

        #calculate merged h value.
        h = float((h1 * Q1 + h2 * Q2) / Q)

        #calculate fluid properties
        prop = multiRP.flsh(u'ph', p1, h, x)

        #calculate Q
        prop[u'Q'] = float(Q)

        #create equipment spec
        equipspec ={u'type': u'FLOWMERGE', u'name': name, u'fluid1': prop_1,
                    u'fluid2': prop_2, u'fluidout': prop}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(u':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop


def pwr_gen(equipspecs):
    u"""Returns the total power generated - the total power consumed

    input:
        equipspecs--equipspec output from equipment in list format"""
    pwr = 0
    for key in equipspecs:
        if u'pwr' in equipspecs[key].keys():
            pwr += equipspecs[key][u'pwr']
    return pwr


def exch_err(equipspecs, criteria=None):
    u"""Returns error margin of exchangers:
    input:
        equipspecs--equipspec output from equipment in list format
        criteria (optional)--acceptance criteria between 0 - 1
    output:
        without criteria
            float no--max. error margin from exchangers
                0 = error free
                1 = 100% error
        with criteria
            True--if error margin <= criteria input
            False--if error margin > criteria input"""
    #raise empty error list
    err = []
    #populate error list
    for key in equipspecs:
        if equipspecs[key][u'type'] == u'EXCHANGER':
            err.append(equipspecs[key][u'errvalue'])
    #Criteria undefined
    if criteria == None:
        #return 0 if equipspecs input does not contain exchanger
        if len(err) == 0:
            return 0
        #return error value
        else: return max(err)
    #Criteria defined through input
    elif 0 <= criteria <= 1:
        #return false if equipspecs input does not contain exchanger
        if len(err) == 0:
            return False
        #return True / False
        else: return max(err) <= criteria
    #criteria not within 0 and 1
    else:
        raise EquipmentinputError(u"input value criteria to be between 0 and " +
                                    u"1 instead of " + unicode(criteria))


def cycletotal(equipspecs, prescript=u'', postscript=u'', mRP=None):
    u"""returns a full overview of total power consumption and all
    equipment details

    input:
        equipspecs--equipspec output from equipment in list format
        prescript--text input before overview
        postscript--text input after overview"""
    #manupulate input equipspecs
    turbines = []
    turbineseps = []
    turbinebleedseps = []
    pumps = []
    compressors = []
    exchangers = []
    subcools = []
    separators = []
    regvalves = []
    flowmerges = []
    for key in equipspecs:
        if equipspecs[key][u'type'] == u'TURBINE':
            turbines.append(equipspecs[key])
        elif equipspecs[key][u'type'] == u'TURBINESEP':
            turbineseps.append(equipspecs[key])
        elif equipspecs[key][u'type'] == u'TURBINEBLEEDSEP':
            turbinebleedseps.append(equipspecs[key])
        elif equipspecs[key][u'type'] == u'PUMP':
            pumps.append(equipspecs[key])
        elif equipspecs[key][u'type'] == u'COMPRESSOR':
            compressors.append(equipspecs[key])
        elif equipspecs[key][u'type'] == u'EXCHANGER':
            exchangers.append(equipspecs[key])
        elif equipspecs[key][u'type'] == u'SUBCOOL':
            subcools.append(equipspecs[key])
        elif equipspecs[key][u'type'] == u'SEPARATOR':
            separators.append(equipspecs[key])
        elif equipspecs[key][u'type'] == u'REGVALVE':
            regvalves.append(equipspecs[key])
        elif equipspecs[key][u'type'] == u'FLOWMERGE':
            flowmerges.append(equipspecs[key])
    #equipment string
    equipstr = u''

    #power string
    powergen = 0
    powercon = 0
    pwrstr = u'*' * 80
    pwrstr += u'{:}{:^80}{:}'.format(u'\n', u'POWER GENERATION', u'\n')

    #turbine
    equipdisplay = True
    for turbine in turbines:
        #power string cont.
        #ensure Turbine is displayed once and only if turbine type in equipspec
        if equipdisplay:
            pwrstr += u'TURBINE\n'
            equipdisplay = False
        pwrstr += _prop(turbine[u'name'], turbine[u'pwr'], u'  Watt')
        powergen += turbine[u'pwr']
        #equipment string cont.
        equipstr += u'*' * 80 + u'\n'
        equipstr += u'{:^80}'.format(u'TURBINE') + u'\n'
        equipstr+= u'{:^80}'.format(turbine[u'name']) + u'\n\n'
        equipstr += _prop(u'Shaft power output: ', turbine[u'pwr'],
                          u'  Watt')
        equipstr += _prop(u'Turbine eff.: ', turbine[u'eff'] * 100, u'%')
        equipstr += _prop(u'Alter. eff.: ', turbine[u'e_eff'] * 100, u'%')
        equipstr += _prop(u'Min. vapor outl.: ', turbine[u'wet'] * 100,
                          u' mol(v)/mol')
        equipstr += _prop(u'add. back press.: ', turbine[u'incr_back_press'],
                          u' kPa')
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID IN', u'\n')
        equipstr += _fluidprop(turbine[u'fluidin'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID OUT', u'\n')
        equipstr += _fluidprop(turbine[u'fluidout'])

    #turbine separator
    equipdisplay = True
    for turbinesep in turbineseps:
        #power string cont.
        #ensure Turbine is displayed once and only if turbine type in equipspec
        if equipdisplay:
            pwrstr += u'TURBINE SEPARATOR\n'
            equipdisplay = False
        pwrstr += _prop(turbinesep[u'name'], turbinesep[u'pwr'], u'  Watt')
        powergen += turbinesep[u'pwr']
        #equipment string cont.
        equipstr += u'*' * 80 + u'\n'
        equipstr += u'{:^80}'.format(u'TURBINE SEPARATOR') + u'\n'
        equipstr+= u'{:^80}'.format(turbinesep[u'name']) + u'\n\n'
        equipstr += _prop(u'Shaft power output: ', turbinesep[u'pwr'],
                          u'  Watt')
        equipstr += _prop(u'Turbine eff.: ', turbinesep[u'eff'] * 100, u'%')
        equipstr += _prop(u'Alter. eff.: ', turbinesep[u'eff'] * 100, u'%')
        equipstr += _prop(u'separator dp: ', turbinesep[u'dp'], u' kPa')
        equipstr += _prop(u'Min. vapor outl.: ', turbinesep[u'wet'] * 100,
                          u' mol(v)/mol')
        equipstr += _prop(u'add. back press.: ', turbinesep[u'incr_back_press'],
                          u' kPa')
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID IN', u'\n')
        equipstr += _fluidprop(turbinesep[u'fluidin'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID OUT', u'\n')
        equipstr += _fluidprop(turbinesep[u'fluidout'])

    #power string cont.
    pwrstr += _prop(u'Subtotal: ', powergen, u'  Watt')
    pwrstr += u'{:}{:^80}{:}'.format(u'\n', u'POWER CONSUMPTION', u'\n')

    #turbine bleed separator
    equipdisplay = True
    for turbinebleedsep in turbinebleedseps:
        #power string cont.
        #ensure Turbine is displayed once and only if turbine type in equipspec
        if equipdisplay:
            pwrstr += u'TURBINE BLEED SEPARATOR\n'
            equipdisplay = False
        pwrstr += _prop(turbinebleedsep[u'name'], turbinebleedsep[u'pwr'],
                        u'  Watt')
        powergen += turbinebleedsep[u'pwr']
        #equipment string cont.
        equipstr += u'*' * 80 + u'\n'
        equipstr += u'{:^80}'.format(u'TURBINE BLEED SEPARATOR') + u'\n'
        equipstr+= u'{:^80}'.format(turbinebleedsep[u'name']) + u'\n\n'
        equipstr += _prop(u'Shaft power output: ', turbinebleedsep[u'pwr'],
                          u'  Watt')
        equipstr += _prop(u'Turbine eff.: ', turbinebleedsep[u'eff'] * 100, u'%')
        equipstr += _prop(u'Alter. eff.: ', turbinebleedsep[u'eff'] * 100, u'%')
        equipstr += _prop(u'separator dp: ', turbinebleedsep[u'dp'], u' kPa')
        equipstr += _prop(u'flow ratio: ',
                          turbinebleedsep[u'Q_ratio'] * 100, u'%')
        equipstr += _prop(u'cond. flow ratio: ',
                          turbinebleedsep[u'Q_ratiobleed'] * 100, u'%')
        equipstr += _prop(u'Min. vapor outl.: ', turbinebleedsep[u'wet'] * 100,
                          u' mol(v)/mol')
        equipstr += _prop(u'add. back press.: ',
                          turbinebleedsep[u'incr_back_press'], u' kPa')
        equipstr += _prop(u'add. bk pr. bleed: ',
                          turbinebleedsep[u'incr_back_press'], u' kPa')
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID IN', u'\n')
        equipstr += _fluidprop(turbinebleedsep[u'fluidin'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID BLEED OUT', u'\n')
        equipstr += _fluidprop(turbinebleedsep[u'fluidbleedout'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID OUT', u'\n')
        equipstr += _fluidprop(turbinebleedsep[u'fluidout'])

    #power string cont.
    pwrstr += _prop(u'Subtotal: ', powergen, u'  Watt')
    pwrstr += u'{:}{:^80}{:}'.format(u'\n', u'POWER CONSUMPTION', u'\n')

    #pump
    equipdisplay = True
    for pump in pumps:
        #power string cont.
        #ensure pump is displayed once and only if pump type in equipspec
        if equipdisplay:
            pwrstr += u'PUMP\n'
            equipdisplay = False
        pwrstr += _prop(pump[u'name'], -pump[u'pwr'], u'  Watt')
        powercon += -pump[u'pwr']
        #equipment string cont.
        equipstr += u'*' * 80 + u'\n'
        equipstr += u'{:^80}'.format(u'PUMP') + u'\n'
        equipstr += u'{:^80}'.format(pump[u'name']) + u'\n\n'
        equipstr += _prop(u'Power consumed: ', -pump[u'pwr'], u'  Watt')
        equipstr += _prop(u'Pump eff.: ', pump[u'eff'] * 100, u'%')
        equipstr += _prop(u'motor eff.: ', pump[u'e_eff'] * 100, u'%')
        equipstr += _prop(u'pump head: ', pump[u'dp'], u' kPa', pump[u'dp'] / 100,
                          u' bar(a)')
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID IN', u'\n')
        equipstr += _fluidprop(pump[u'fluidin'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID OUT', u'\n')
        equipstr += _fluidprop(pump[u'fluidout'])

    #compressor
    equipdisplay = True
    for compressor in compressors:
        #power string cont.
        #ensure pump is displayed once and only if pump type in equipspec
        if equipdisplay:
            pwrstr += u'COMPRESSOR\n'
            equipdisplay = False
        pwrstr += _prop(compressor[u'name'], -compressor[u'pwr'], u'  Watt')
        powercon += -compressor[u'pwr']
        #equipment string cont.
        equipstr += u'*' * 80 + u'\n'
        equipstr += u'{:^80}'.format(u'COMPRESSOR') + u'\n'
        equipstr += u'{:^80}'.format(compressor[u'name']) + u'\n\n'
        equipstr += _prop(u'Power consumed: ', -compressor[u'pwr'], u'  Watt')
        equipstr += _prop(u'Compr. eff.: ', compressor[u'eff'] * 100, u'%')
        equipstr += _prop(u'Motor eff.: ', compressor[u'e_eff'] * 100, u'%')
        equipstr += _prop(u'compressor head: ', compressor[u'dp'], u' kPa')
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID IN', u'\n')
        equipstr += _fluidprop(compressor[u'fluidin'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID OUT', u'\n')
        equipstr += _fluidprop(compressor[u'fluidout'])

    #power string cont.
    pwrstr += _prop(u'Subtotal: ', powercon, u'  Watt')
    pwrnett = powergen - powercon
    if pwrnett == scipy.inf:
        pwrnett = -scipy.inf
    pwrstr += u'{:}{:^80}{:}'.format(u'\n', u'NETT POWER GENERATION', u'\n')
    pwrstr += _prop(u'Nett power gen.: ', pwrnett, u'  Watt')

    #exchanger
    for exchanger in exchangers:
        #equipment string cont.
        equipstr += u'*' * 80 + u'\n'
        equipstr += u'{:^80}'.format(u'EXCHANGER') + u'\n'
        equipstr += u'{:^80}'.format(exchanger[u'name']) + u'\n\n'
        equipstr += _prop(u'Cooling duty (cold): ', exchanger[u'E_cold'],
                          u'  Watt')
        equipstr +=_prop(u'Cooling duty (hot): ', exchanger[u'E_hot'],
                         u'  Watt')
        equipstr += _prop(u'Error value: ', exchanger[u'errvalue'] * 100,
                          u'%')
        equipstr += _prop(u'Delta temp.: ', exchanger[u'dt'], u'  Kelvin')
        #if fluidin = cold
        if exchanger[u'fluidin'][u't'] < exchanger[u'fluidin_contra'][u't']:
            equipstr += _prop(u'Delta press. cold: ', exchanger[u'dp'], u'  kPa')
            equipstr += _prop(u'Delta press. hot: ', exchanger[u'dp_contra'],
                              u'  kPa')
            equipstr += u'{:}{:^80}{:}'.format(u'\n', u'COLD FLUID IN', u'\n')
            equipstr += _fluidprop(exchanger[u'fluidin'])
            equipstr += u'{:}{:^80}{:}'.format(u'\n', u'COLD FLUID OUT', u'\n')
            equipstr += _fluidprop(exchanger[u'fluidout'])
            equipstr += u'{:}{:^80}{:}'.format(u'\n', u'HOT FLUID IN', u'\n')
            equipstr += _fluidprop(exchanger[u'fluidin_contra'])
            equipstr += u'{:}{:^80}{:}'.format(u'\n', u'HOT FLUID OUT', u'\n')
            equipstr += _fluidprop(exchanger[u'fluidout_contra'])
        #if fluidin = hot
        else:
            equipstr += _prop(u'Delta press. cold: ', exchanger[u'dp_contra'],
                              u'  kPa')
            equipstr += _prop(u'Delta press. hot: ', exchanger[u'dp'], u'  kPa')
            equipstr += u'{:}{:^80}{:}'.format(u'\n', u'COLD FLUID IN', u'\n')
            equipstr += _fluidprop(exchanger[u'fluidin_contra'])
            equipstr += u'{:}{:^80}{:}'.format(u'\n', u'COLD FLUID OUT', u'\n')
            equipstr += _fluidprop(exchanger[u'fluidout_contra'])
            equipstr += u'{:}{:^80}{:}'.format(u'\n', u'HOT FLUID IN', u'\n')
            equipstr += _fluidprop(exchanger[u'fluidin'])
            equipstr += u'{:}{:^80}{:}'.format(u'\n', u'HOT FLUID OUT', u'\n')
            equipstr += _fluidprop(exchanger[u'fluidout'])
        #add graph
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLOW - TEMP.', u'\n')
        equipstr += u'{:^80}{:}'.format(u'GRAPHIC', u'\n')
        try:
            graph = _exchgraph(exchanger[u'fluidin'],
                               exchanger[u'fluidout'],
                               exchanger[u'fluidin_contra'],
                               exchanger[u'fluidout_contra'], mRP=mRP)
        except:
            graph = u'error in generating graphics at def cycletotal\n'
        equipstr += graph

    #subcool
    for subcool in subcools:
        #equipment string cont.
        equipstr += u'*' * 80 + u'\n'
        equipstr += u'{:^80}'.format(u'SUBCOOL') + u'\n'
        equipstr += u'{:^80}'.format(subcool[u'name']) + u'\n\n'
        equipstr += _prop(u'Delta temp.: ', subcool[u'dt'], u'  Kelvin')
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID IN', u'\n')
        equipstr += _fluidprop(subcool[u'fluidin'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID OUT', u'\n')
        equipstr += _fluidprop(subcool[u'fluidout'])

    #separator
    for separator in separators:
        #equipment strin cont.
        equipstr += u'*' * 80 + u'\n'
        equipstr += u'{:^80}'.format(u'SEPARATOR') + u'\n'
        equipstr += u'{:^80}'.format(separator[u'name']) + u'\n\n'
        equipstr += _prop(u'Delta press.: ', separator[u'dp'], u'  kPa')
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID IN', u'\n')
        equipstr += _fluidprop(separator[u'fluidin'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'VAPOR OUT', u'\n')
        equipstr += _fluidprop(separator[u'fluidvap'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'LIQUID OUT', u'\n')
        equipstr += _fluidprop(separator[u'fluidliq'])

    #regulating valve
    for regvalve in regvalves:
        #equipment strin cont.
        dp = regvalve[u'fluidin'][u'p'] - regvalve[u'p_out']
        equipstr += u'*' * 80 + u'\n'
        equipstr += u'{:^80}'.format(u'REGULATING VALVE') + u'\n'
        equipstr += u'{:^80}'.format(regvalve[u'name']) + u'\n\n'
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID IN', u'\n')
        equipstr += _fluidprop(regvalve[u'fluidin'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID OUT', u'\n')
        equipstr += _fluidprop(regvalve[u'fluidout'])

    #2-flow merge
    for flowmerge in flowmerges:
        #equipment strin cont.
        equipstr += u'*' * 80 + u'\n'
        equipstr += u'{:^80}'.format(u'FLOWMERGE') + u'\n'
        equipstr += u'{:^80}'.format(flowmerge[u'name']) + u'\n\n'
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID 1 IN', u'\n')
        equipstr += _fluidprop(flowmerge[u'fluid1'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID 2 IN', u'\n')
        equipstr += _fluidprop(flowmerge[u'fluid2'])
        equipstr += u'{:}{:^80}{:}'.format(u'\n', u'FLUID OUT', u'\n')
        equipstr += _fluidprop(flowmerge[u'fluidout'])

    #last line
    equipstr += u'*' * 80 + u'\n'

    #return total string
    return prescript + u'\n' + pwrstr + equipstr + postscript


#for testing only
def _test(prop, mRP=None):

    equipdct = {}

    print u'nett power before equipment generation'
    print pwr_gen(equipdct)

    print u'run subcool'
    prop, equipdct[u'subcool_01'] = subcool(prop, 4.5, u'subcool-01')

    print u'run pump'
    prop03, equipdct[u'pump_01'] = pump(prop, 1000, 0.85, 0.94, name=u'pump-01')

    #exch_err(0.4)
    print u'run exchanger simple'
    prop04, equipdct[u'exchanger_01'] = exchanger(prop03, 25, 4,
                                                 name=u'exchanger-01')

    print u'run separator'
    prop04[u't'] += 100
    prop01 = multiRP.flsh(u'tp', prop04[u't'], prop04[u'p'], prop04[u'x'])
    prop01[u'Q'] = 1
    propliq, prop02, equipdct[u'separator_01'] = separator(prop01, 50,
                                                          u'separator-01')

    print u'run reg_valve'
    prop_reg_valve, equipdct[u'regvalve_01'] = reg_valve(propliq, 100,
                                                        u'regvalve-01')

    print u'run turbine'
    prop, equipdct[u'turbine_01'] = turbine(prop02, 100, 0.85, 0.96,
                                           wet=0.97, name=u'Turbine-01')

    print u'run flowmerge'
    prop_fm, equipdct[u'flowmerge_01'] = flowmerge(prop, prop_reg_valve,
                                                  u'flowmerge-01')

    print u'run turbine_separator'
    prop, equipdct[u'turbine_separator_01'] = turbine_separator(prop01,
        100, 50, 0.85, wet=0.95, name=u'turbine_separator-01')

    print u'run compressor & exchanger inletflow'
    prop05, equipdct[u'compressor_01']= compressor(prop02, 3000, 0.8,
                                                 0.94, name=u'compressor-01')

    print u'run exchanger'
    prop06, equipdct[u'exchanger_01'] = exchanger(prop05, 25, 4, prop03, 25,
                                                 prop04, name=u'Exchanger-01')

    #exch_err(0.1)
    print u'run exchanger otherflow outlet'
    prop07, equipdct[u'exchanger_01'] = exchanger(prop03, 25, 4, prop05, 25,
                                                 prop06, name=u'Exchanger-01')

    if abs((((prop05[u'h'] - prop06[u'h']) * prop05[u'Q']) -
             ((prop07[u'h'] - prop03[u'h']) * prop03[u'Q'])) /
            ((prop05[u'h'] - prop06[u'h']) * prop05[u'Q'])) < 0.00000001:
        print u'exchanger tolerance test acceptable\n'
    else:
        print u'exchanger tolerance test failed\n'

    print u'exchanger err value'
    print exch_err(equipdct)
    print exch_err(equipdct, 0.02)


    print u'nett power'
    print pwr_gen(equipdct)

    return cycletotal(equipdct, prescript=u'start', postscript=u'end')

if __name__ == u'__main__':
    #add module file path to python sys path
    import equipment as _filename
    _filename = (os.path.dirname(_filename.__file__))
    sys.path.append(_filename)

    #examples and test setup
    #setup fluid
    prop = multiRP.setup(u'def', u'water', u'ammonia')
    prop[u'x'] = [0.7, 0.3]
    prop[u't'] = 273.15 + 35
    prop[u'Q'] = 1

    #single processing
    SetMultiprocessing.off()
    print _test(prop)

    #multiprocessing
    SetMultiprocessing.on()
    print _test(prop)

    print u'test'
    print SetMultiprocessing()

