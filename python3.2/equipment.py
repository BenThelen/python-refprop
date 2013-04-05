#-------------------------------------------------------------------------------
#Name:                  equipment
#Purpose:               Calculate equipment fluid output properties
#                       using refprop module
#
#Author:                Thelen, B.J.
#                       thelen_ben@yahoo.com
#-------------------------------------------------------------------------------

'''This module usage the refprop module to calculate fluid output
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
import os, sys, scipy, math, decimal, copy
import multiRP
import multiprocessing as mp

#classes
class EquipmentError(Exception):
    'General EquipmentError for python module'
    pass

class EquipmentinputError(EquipmentError):
    'Equipment input Error'
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class EquipmentinfeasibleError(EquipmentError):
    'Equipment input Error'
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class SetMultiprocessing:
    'Return _setmultiprocessing status (on / off)'
    def __repr__(self):
        if not '_setmultiprocessing' in globals(): SetMultiprocessing.on()
        return _setmultiprocessing
    @staticmethod
    def on():
        'Sets SetMultiprocessing on, allows equipment multiprocessing'
        global _setmultiprocessing
        _setmultiprocessing = 'on'
        return _setmultiprocessing
    @staticmethod
    def off():
        'Sets SetMultiprocessing off, single processing for equipment module'
        global _setmultiprocessing
        _setmultiprocessing = 'off'
        return _setmultiprocessing


#suporting functions
def _inputerrorcheck(deflocals, equipment):
    '''check valid def input or else raise error'''
    for key in deflocals.keys():
        #check input prop
        if key[:4] == 'prop':
            if deflocals[key] != None:
                if not deflocals[key].__class__ == dict:
                    raise EquipmentinputError("expect dict input for 'prop' " +
                                               "instead of " +
                                               str(deflocals['prop'].__class__))

            #check input prop['t']
            if ['subcool', 'exchanger'].__contains__(equipment):
                if deflocals[key] != None:
                    if not 't' in deflocals[key]:
                        raise EquipmentinputError('prop dict key "t" ' +
                                                   'required as input')

            #check input prop['h']
            if ['separator', 'reg_valve', 'turbine', 'flowmerge',
                'turbine_separator',
                'turbine_bleed_sep'].__contains__(equipment):
                if deflocals[key] != None:
                    if not 'h' in deflocals[key]:
                        raise EquipmentinputError('prop dict key "h" ' +
                                                   'required as input')

            #check input prop['x']
            if ['subcool', 'pump', 'separator', 'reg_valve', 'turbine',
                'flowmerge', 'turbine_separator', 'compressor',
                'turbine_bleed_sep', 'exchanger'].__contains__(equipment):
                if deflocals[key] != None:
                    if not 'x' in deflocals[key]:
                        raise EquipmentinputError('prop dict key "x" ' +
                                                   'required as input')

            #check input prop['Q']
            if ['pump', 'separator', 'turbine', 'flowmerge',
            'turbine_separator', 'compressor', 'turbine_bleed_sep',
            'exchanger'].__contains__(equipment):
                if deflocals[key] != None:
                    if not 'Q' in deflocals[key]:
                        raise EquipmentinputError('prop dict key "Q" ' +
                                                   'required as input')
                    elif deflocals[key]['Q'] < 0:
                        raise EquipmentinputError('"Q" value to be > 0')

            #check input prop['s']
            if ['pump', 'turbine', 'compressor'].__contains__(equipment):
                if deflocals[key] != None:
                    if not 's' in deflocals[key]:
                        raise EquipmentinputError('prop dict key "s" ' +
                                                   'required as input')

            #check input prop['p']
            if ['pump', 'separator','flowmerge', 'turbine_separator',
                'compressor', 'exchanger',
                'turbine_bleed_sep'].__contains__(equipment):
                if deflocals[key] != None:
                    if not 'p' in deflocals[key]:
                        raise EquipmentinputError('prop dict key "p" ' +
                                                   'required as input')

            #check equal prop_1 and prop_2['p']
            if ['flowmerge'].__contains__(equipment):
                if 'prop_2' in deflocals.keys() \
                and deflocals['prop_1']['p'] != deflocals['prop_2']['p']:
                    raise EquipmentinputError('pressure value to be equal ' +
                                               'for prop_1 and prop_2')

            #check equal prop_1 and prop_2['hfld']
            if ['flowmerge'].__contains__(equipment):
                if 'prop_2' in deflocals.keys() \
                and deflocals['prop_1']['hfld'] != deflocals['prop_2']['hfld']:
                    raise EquipmentinputError('hfld value to be equal for ' +
                                               'prop_1 and prop_2')

        #check input dt
        if (['subcool'].__contains__(equipment) \
        or ['exchanger'].__contains__(equipment)) \
        and key == 'dt':
            if not deflocals[key].__class__ == float \
            and not deflocals[key].__class__ == int:
                raise EquipmentinputError("expect float/int input for 'dt' " +
                                           'instead of ' +
                                           str(deflocals[key].__class__))
            if deflocals['dt'] < 0:
                raise EquipmentinputError("'dt' value shall be positive")

        #check input eff
        if ['pump', 'turbine', 'turbine_separator', 'turbine_bleed_sep',
            'compressor'].__contains__(equipment) \
            and (key == 'eff' \
                 or key == 'e_eff'):
            if not deflocals[key].__class__ == float \
            and not 0 <= key <= 1:
                raise EquipmentinputError("expect float value between 0 and 1 "
                                           "for '" + key + "' instead of " +
                                           str(deflocals[key]))

        #check input Q_ratio
        if ['turbine_bleed_sep'].__contains__(equipment) \
           and (key == 'Q_ratio' or key == 'Q_ratiobleed'):
            if not deflocals[key].__class__ == float \
            and not 0 < key < 1:
                raise EquipmentinputError("expect float value between 0 and 1 "
                                           "for '" + key + "' instead of " +
                                           str(deflocals[key]))

        #check input dp
        if ['pump', 'separator', 'turbine_separator', 'compressor', 'exchanger',
            'turbine_bleed_sep'].__contains__(equipment) and key == 'dp':
            if not deflocals[key].__class__ == float \
            and not deflocals[key].__class__ == int:
                raise EquipmentinputError("expect float/int input for 'dp' " +
                                           'instead of ' +
                                           str(deflocals[key].__class__))
            if deflocals['dp'] < 0:
                raise EquipmentinputError("'dp' value shall be positive")

        #check input dp_contra
        if ['exchanger'].__contains__(equipment) and key == 'dp_contra':
            if not deflocals[key].__class__ == float \
            and not deflocals[key].__class__ == int:
                raise EquipmentinputError("expect float/int input for '" +
                                           "'dp_contra' instead of " +
                                           str(deflocals[key].__class__))
            if deflocals['dp_contra'] < 0:
                raise EquipmentinputError("'dp_contra' value shall be positive")

        #check input p_out
        if ['reg_valve', 'turbine', 'turbine_bleed_sep',
             'turbine_separator'].__contains__(equipment) \
             and (key == 'p_out' or key == 'p_outbleed'):
            if not deflocals[key].__class__ == float \
            and not deflocals[key].__class__ == int:
                raise EquipmentinputError("expect float/int input for '" +
                                           key + "' instead of " +
                                           str(deflocals[key].__class__))


def _prop(name, value, unit, valueSI=None, unitSI=None):
    """print out fluid property in organised format

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
        fldob = str('{0:22}' + '{1:>13,.3f}' + '{2:18}\n').format(name, value,
                                                                  '')
    if valueSI == None:
        #format print out
        fldob = str('{0:22}' + '{1:>13,.3f}' + '{2:18}\n').format(name, value,
                                                                  unit)
    elif value == scipy.inf or value == -scipy.inf:
        #format print out
        fldob = str('{0:22}' + '{1:>13,.3f}' + '{2:19}' + '{3:>13,.3f}' +
                    '{4}\n').format(name, value, '', valueSI, '')
    else:
        #format print out
        fldob = str('{0:22}' + '{1:>13,.3f}' + '{2:18}' + '{3:>13,.3f}' +
                    '{4}\n').format(name, value, unit, valueSI, unitSI)
    return fldob


def _fluidprop(prop_in):
    """print all fluid properties in organised format

    input:
        prop_in--fluid characteristics to be plotted
    output:
        prpstr--fluid property in organized string format"""
    prpstr = ''
    if 'hfld' in prop_in.keys():
        multiRP.resetup(prop_in)
        if 'x' in prop_in.keys():
            x = []
            for each in range(len(prop_in['x'])):
                x.insert(each, float(prop_in['x'][each]))
        if 'xliq' in prop_in.keys():
            xliq = []
            for each in range(len(prop_in['xliq'])):
                xliq.insert(each, float(prop_in['xliq'][each]))
        if 'xvap' in prop_in.keys():
            xvap = []
            for each in range(len(prop_in['xvap'])):
                xvap.insert(each, float(prop_in['xvap'][each]))
        #determine mole weigth composition
        #assigned [0] by echanger to correct for wmix = 0
        if x == [0]:
            wmix = 0
        else:
            wmix = multiRP.xmass(x)['wmix']
        if 'q' in prop_in.keys():
            if 0 < prop_in['q'] < 1:
                if 'xliq' in prop_in.keys():
                    wmixliq = multiRP.xmass(xliq)['wmix']
                if 'xvap' in prop_in.keys():
                    wmixvap = multiRP.xmass(xvap)['wmix']
        #fluid components display
        #correction on exchanger assigned phantom fluid
        if 'hfld' in prop_in.keys():
            for each in range(len(prop_in['hfld'])):
                #obtain mole weigth individual comp.
                wmm = multiRP.info(each + 1)['wmm']
                prpstr += _prop(str(prop_in['hfld'][each]),
                                x[each] * 100,
                                '% mol(i)/mol',
                                x[each] / wmix * wmm * 100,
                                '% g(i)/g')
        multiRP.resetup(prop_in)
        #liquid components display
        if 'q' in prop_in.keys():
            if 0 < prop_in['q'] < 1:
                for each in range(len(prop_in['hfld'])):
                    #obtain mole weigth individual comp.
                    wmm = multiRP.info(each + 1)['wmm']
                    prpstr += _prop(str(prop_in['hfld'][each]) + ' (liq)',
                                    xliq[each] * 100,
                                    '% mol(i)/mol',
                                    xliq[each] / wmixliq * wmm * 100,
                                    '% g(i)/g')
        multiRP.resetup(prop_in)
        #vapor components display
        if 'q' in prop_in.keys():
            if 0 < prop_in['q'] < 1:
                for each in range(len(prop_in['hfld'])):
                    #obtain mole weigth individual comp.
                    wmm = multiRP.info(each + 1)['wmm']
                    prpstr += _prop(str(prop_in['hfld'][each]) + ' (vap)',
                                    xvap[each] * 100,
                                    '% mol(i)/mol',
                                    xvap[each] / wmixvap * wmm * 100,
                                    '% g(i)/g')
        multiRP.resetup(prop_in)
        #Flow display
        if 'Q' in prop_in.keys():
            prpstr += _prop('Flow: ',
                            prop_in['Q'],
                            '  mol/s',
                            prop_in['Q'] * wmix * 3600 / 1000,
                            '  kg/h')
        #mol weight display
        #correction for exchanger phantom fluid
        if not wmix == 0:
            prpstr += _prop('molecular weight: ',
                            wmix,
                            '  g/mol')
        #pressure display
        if 'p' in prop_in.keys():
            prpstr += _prop('pressure: ',
                            prop_in['p'],
                            '  kPa(a)',
                            prop_in['p'] / 100,
                            '  bar(a)')
        #temp display
        if 't' in prop_in.keys():
            prpstr += _prop('temperature: ',
                            prop_in['t'],
                            '  Kelvin',
                            prop_in['t'] - 273.15,
                            '  Celsius')
        #Density display
        if 'D' in prop_in.keys():
            prpstr += _prop('Density: ',
                            prop_in['D'],
                            '  mol/L',
                            prop_in['D'] * wmix,
                            '  kg/m3')
        #liq Density display
        if 'q' in prop_in.keys():
            if 0 < prop_in['q'] < 1:
                if 'Dliq' in prop_in.keys():
                    prpstr += _prop('Density (liq): ',
                                    prop_in['Dliq'],
                                    '  mol/L',
                                    prop_in['Dliq'] * wmixliq,
                                    '  kg/m3')
        #vap density display
        if 'q' in prop_in.keys():
            if 0 < prop_in['q'] < 1:
                if 'Dvap' in prop_in.keys():
                    prpstr += _prop('Density (vap): ',
                                    prop_in['Dvap'],
                                    '  mol/L',
                                    prop_in['Dvap'] * wmixvap,
                                    '  kg/m3')
        #enthalpy display
        if 'h' in prop_in.keys() and prop_in['h'] != 0:
            prpstr += _prop('Enthalpy: ',
                            prop_in['h'],
                            '  J/mol',
                            prop_in['h'] / wmix,
                            '  kJ/kg')
        #entropy display
        if 's' in prop_in.keys() and prop_in['h'] != 0:
            prpstr += _prop('Entropy: ',
                            prop_in['s'],
                            '  J/(molK)',
                            prop_in['s'] / wmix,
                            '  kJ/kgK')
        #phase display
        try:
            phase = multiRP.getphase(prop_in)
            #quality display
            if phase == "2 phase":
                prpstr += _prop('vapor quality: ', prop_in['q'], '  mol(v)/mol',
                                prop_in['q'] * wmixvap / wmix, '  g(v)/g')
            #liquid (subcool) display
            elif phase == "liquid":
                t = multiRP.satp(prop_in['p'], x, 1)['t']
                prpstr += _prop('subcool: ', t - prop_in['t'], '  Kelvin',
                                t - prop_in['t'], '  Celsius')
            #saturated liquid display
            elif phase == "saturated liquid":
                prpstr += 'saturated liquid\n'
            #saturated vapor display
            elif phase == "saturated vapor":
                prpstr += 'saturated vapor\n'
            #superheated vapor display
            elif phase == "vapor":
                t = multiRP.satp(prop_in['p'], prop_in['x'], 2)['t']
                prpstr += _prop('superheated vapor: ', prop_in['t'] - t, '  Kelvin',
                                prop_in['t'] - t, '  Celsius')
            #superheated gas display
            elif phase == "gas":
                t = multiRP.satp(prop_in['p'], prop_in['x'], 2)['t']
                prpstr += _prop('superheated gas: ', prop_in['t'] - t, '  Kelvin',
                                prop_in['t'] - t, '  Celsius')
            #compressible liquid display
            elif phase == "compressible liquid":
                pcrit = multiRP.critp(prop_in['x'])['pcrit']
                prpstr += _prop('compres. lqd: ', prop_in['p'] - pcrit,
                                '  kPa(a)', (prop_in['p'] - pcrit) / 100,
                                '  bar(a)')
            #supercritical fluid display
            elif phase == "Supercritical fluid":
                prop = multiRP.critp(prop_in['x'])
                pcrit = prop['pcrit']
                tcrit = prop['tcrit']
                prpstr += _prop('Supercritical: ', prop_in['p'] - pcrit,
                                '  kPa(a)', (prop_in['p'] - pcrit) / 100,
                                '  bar(a)')
                prpstr += _prop('', prop_in['t'] - tcrit, '  Kelvin',
                                prop_in['t'] - tcrit, '  Celsius')
        except:
            pass
    #return fluid object
    return prpstr


def _exchgraph(prop_in, prop_out, contra_prop_in, contra_prop_out, mRP=None):
    """return a string representing the exchanger in /out fluids temp graph"""
    try:
        lines = 30
        position = 80
        c, h = [], []
        #define cold and hot fluids
        if prop_in['t'] < contra_prop_in['t']:
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
        h_cold_step = (cold_out['h'] - cold_in['h'])/(position-1)
        p_cold_step = (cold_out['p'] - cold_in['p'])/(position-1)
        h_hot_step = (hot_in['h'] - hot_out['h'])/(position-1)
        p_hot_step = (hot_in['p'] - hot_out['p'])/(position-1)

        #check if multiprocessing is allowed
        if SetMultiprocessing().__repr__() == 'on':
            ##start with multiprocessing##
            if mRP == None:
                #create if not exist
                mRP = multiRP.multirefprop()#initiate multirefprop

            ##create children
            for pos in range(position):
                c.append(mRP['process'](target=_coldhot,
                                        args=(cold_in['p'] + p_cold_step * pos,
                                              cold_in['h'] + h_cold_step * pos,
                                              cold_in['x']),
                                        kwargs={'prop':cold_in, 'mRP':mRP}))
                h.append(mRP['process'](target=_coldhot,
                                        args=(hot_out['p'] + p_hot_step * pos,
                                              hot_out['h'] + h_hot_step * pos,
                                              hot_in['x']),
                                        kwargs={'prop':hot_in, 'mRP':mRP}))
            fldlist = c + h

            #run the multiprocessing list
            multiRP.run_mRP(fldlist)
            #obtain results from multiprocessing list
            templist = []
            for each in fldlist:
                if each.name in mRP['result']:
                    if mRP['result'][each.name].__class__ == dict \
                    and 't' in mRP['result'][each.name]:
                        templist.append(mRP['result'][each.name]['t'])
                    else:
                        templist.append(mRP['result'][each.name])
            ##end with multiprocessing##

        #check for single processing
        elif SetMultiprocessing().__repr__() == 'off':
            #loop to get each position
            for pos in range(position):
                c.append(_coldhot(cold_in['p'] + p_cold_step * pos,
                                  cold_in['h'] + h_cold_step * pos,
                                  cold_in['x'],
                                  prop=cold_in))
                h.append(_coldhot(hot_out['p'] + p_hot_step * pos,
                                  hot_out['h'] + h_hot_step * pos,
                                  hot_in['x'],
                                  prop=hot_in))
            #create complete temp. list for each position
            templist = [each if each == None else each['t'] for each in (c + h)]

        #round of due to graph display hickups
        for each in range(len(templist)):
            if templist[each] != None:
                templist[each] = round(templist[each], 10)

        #split templist (hot and cold)
        c = templist[:int(len(templist)/2)]
        h = templist[int(len(templist)/2):]

        #determine temperature range and split evenly for 30 lines
        dt_step = (hot_in['t'] - cold_in['t']) / (lines - 1)
        #create lines
        graph = ''
        for line in range(lines):
            #create positions within line
            for pos in range(position):
                #if both h and c then plot X
                if c[pos] != None and h[pos] != None \
                and hot_in['t'] - dt_step * (line + 1) < \
                c[pos] <= \
                hot_in['t'] - dt_step * line \
                and hot_in['t'] - dt_step * (line + 1) < \
                h[pos] <= \
                hot_in['t'] - dt_step * line:
                    graph += 'X'
                #if c only then plot C
                elif c[pos] != None \
                and hot_in['t'] - dt_step * (line + 1) < \
                c[pos] <= \
                hot_in['t'] - dt_step * line:
                    graph += 'C'
                #if h only then plot H
                elif h[pos] != None \
                and hot_in['t'] - dt_step * (line + 1) < \
                h[pos] <= \
                hot_in['t'] - dt_step * line:
                    graph += 'H'
                #else plot space
                else: graph += ' '
            #next line
            graph += '\n'
        return graph
    except:
        return 'error in generating graphics at def _exchgraph\n'

def _coldhot(p, h, x, prop=None, mRP=None):
    '''cold / hot function for _exchgraph'''
    try:
        return multiRP.flsh('ph', p, h, x, prop=prop, mRP=mRP)
    except multiRP.RefpropError:
        #check for multiprocessing
        if SetMultiprocessing().__repr__() == 'on':
            multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = None
        #check for single processing
        elif SetMultiprocessing().__repr__() == 'off':
            return None



#primary functions
def pump(prop_in, dp, eff=1.0, e_eff=1.0, name=False, mRP=None):
    """Calculate fluid properties at pump discharge.
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
    _inputerrorcheck(locals(), 'pump')
    multiRP.resetup(prop_in, mRP=mRP)
    p = prop_in['p'] + dp
    s = prop_in['s']
    x = prop_in['x']
    flshcalc = False
    if 'h' in prop_in: h = prop_in['h']
    else:
        #inhibit error reporting
        if str(multiRP.SetErrorDebug()) == 'on':
            sed = multiRP.SetErrorDebug.on
            multiRP.SetErrorDebug.off()
        elif str(multiRP.SetErrorDebug()) == 'off':
            sed = multiRP.SetErrorDebug.off
        try: #use flsh1 calculation for speed increase
            h = multiRP.psliq(p, s, x)['h']
            prop.update(sed())
        except multiRP.RefpropError:
            sed()
            h = multiRP.flsh('ps', p, s, x)['h']
            flshcalc = True
    Q = prop_in['Q']

    #calculate enthalpy at 100% efficiency
    #inhibit error reporting
    if str(multiRP.SetErrorDebug()) == 'on':
        sed = multiRP.SetErrorDebug.on
        multiRP.SetErrorDebug.off()
    elif str(multiRP.SetErrorDebug()) == 'off':
        sed = multiRP.SetErrorDebug.off
    try: #use flsh1 calculation for speed increase
        if flshcalc: raise multiRP.RefpropError
        prop = multiRP.psliq(p, s, x)
        prop.update(sed())
    except multiRP.RefpropError:
        sed()
        prop = multiRP.flsh('ps', p, s, x)
        flshcalc = True

    #correct h for efficiency losses
    h += (prop['h'] - h) / eff

    #calculate properties with corrected h and discharge p
    #inhibit error reporting
    if str(multiRP.SetErrorDebug()) == 'on':
        sed = multiRP.SetErrorDebug.on
        multiRP.SetErrorDebug.off()
    elif str(multiRP.SetErrorDebug()) == 'off':
        sed = multiRP.SetErrorDebug.off
    try: #use flsh1 calculation for speed increase
        if flshcalc: raise multiRP.RefpropError
        prop = multiRP.phliq(p, h, x)
        prop.update(sed())
    except multiRP.RefpropError:
        sed()
        prop = multiRP.flsh('ph', p, h, x)

    #restore Q value
    prop['Q'] = Q

    #power consumption
    pwr = ((prop['h'] - prop_in['h']) * Q) / e_eff

    equipspec ={'type': 'PUMP', 'name': name, 'dp': dp, 'eff': eff,
                'fluidin': prop_in, 'fluidout': prop, 'pwr': -pwr,
                'e_eff':e_eff}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop


def exchanger(prop_in, dp, dt, prop_contra_in={}, dp_contra=0,
              prop_contra_out={}, name=False, mRP=None):
    """Calculate fluid out properties of exchanger.
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
    if (not 't' in locals()['prop_contra_out']
         or not 'Q' in locals()['prop_contra_out']
         or prop_contra_out['Q'] <= 0):
        prop_contra_out['t'] = 0
    if not 'Q' in locals()['prop_contra_in']:
        prop_contra_out['Q'] = 0
        prop_contra_in['Q'] = 0
    else: prop_contra_out['Q'] = prop_contra_in['Q']
    if not 'p' in locals()['prop_contra_in']:
        prop_contra_out['p'] = 0
        prop_contra_in['p'] = 0
    else: prop_contra_out['p'] = prop_contra_in['p'] - dp_contra
    if not 'x' in locals()['prop_contra_in']:
        prop_contra_out['x'] = [0]
        prop_contra_in['x'] = [0]
    else: prop_contra_out['x'] = prop_contra_in['x']
    if not 't' in locals()['prop_contra_in']:
        prop_contra_in['t'] = 0
    if not 'h' in locals()['prop_contra_in']:
        prop_contra_in['h'] = 0
    if not 'h' in locals()['prop_contra_out']:
        prop_contra_out['h'] = 0

    _inputerrorcheck(locals(), 'exchanger')

    #define variables with values
    Q_in = prop_in['Q']
    p_in = prop_in['p']
    t_in = prop_in['t']
    x_in = prop_in['x']
    if not 'h' in prop_in:
        multiRP.resetup(prop_in, mRP=mRP)
        ########################################################################
        ############ actual h should be present for pure fluids ################
        ########################################################################
        h_in = multiRP.flsh('tp', t_in, p_in, x_in)['h']
    else: h_in = prop_in['h']
    p_out = p_in - dp
    if p_out < 0: p_out = 0

    #define 0 contra flow 'for initiating of calc'
    if prop_contra_out['t'] == 0 or abs(t_in - prop_contra_in['t']) < dt:
        multiRP.resetup(prop_in, mRP=mRP)
        prop = multiRP.flsh('ph', p_out, h_in, x_in)

    #define Dh based on calculations
    else:
        #define cold and hot prop stream
        if prop_in['t'] <= prop_contra_in['t']:
            prop_cold_in = copy.copy(prop_in)
            prop_hot_in = copy.copy(prop_contra_in)
            dp_cold = dp
            dp_hot = dp_contra
        elif prop_in['t'] > prop_contra_in['t']:
            prop_hot_in = copy.copy(prop_in)
            prop_cold_in = copy.copy(prop_contra_in)
            dp_cold = dp_contra
            dp_hot = dp
        Q_cold_in = prop_cold_in['Q']
        t_cold_in = prop_cold_in['t']
        p_cold_in = prop_cold_in['p']
        x_cold_in = prop_cold_in['x']
        if not 'q' in prop_cold_in: q_cold_in = -1
        else: q_cold_in = prop_cold_in['q']
        if not 'h' in prop_cold_in:
            multiRP.resetup(prop_cold_in, mRP=mRP)
        else: h_cold_in = prop_cold_in['h']
        Q_hot_in = prop_hot_in['Q']
        t_hot_in = prop_hot_in['t']
        p_hot_in = prop_hot_in['p']
        x_hot_in = prop_hot_in['x']
        if not 'q' in prop_hot_in: q_hot_in = 2
        else: q_hot_in = prop_hot_in['q']
        if not 'h' in prop_hot_in:
            multiRP.resetup(prop_hot_in, mRP=mRP)
        else: h_hot_in = prop_hot_in['h']

        #define max / min temp of fluid allowed for refprop
        #cold fluid (t_max)
        multiRP.resetup(prop_cold_in, mRP=mRP)
        if 'setmod' in prop_cold_in and 'htype' in prop_cold_in['setmod']:
            tmax_cold = (multiRP.limits(x_cold_in,
                prop_cold_in['setmod']['htype'])['tmax'] - dt)
        else: tmax_cold = multiRP.limits(x_cold_in)['tmax'] - dt

        #hot fluid (t_min)
        multiRP.resetup(prop_hot_in)
        if 'setmod' in prop_hot_in and 'htype' in prop_hot_in['setmod']:
            #check freezing point
            #IMPROVEMENT OPTION, refprop is incapable to calc. freezing points
            #at various mixtures, also the FLASH calculations are limited
            tmin_ht1 = max([multiRP.limitk(htype=prop_hot_in['setmod']['htype'],
                                            icomp=(i+1))['tmin'] \
                            for i in range(prop_hot_in['nc'])])
            #check limits
            tmin_ht2 = (multiRP.limits(x_hot_in,
                                       prop_hot_in['setmod']['htype'])['tmin'])
        else:
            #check freezing point
            #IMPROVEMENT OPTION, refprop is incapable to calc. freezing points
            #at various mixtures, also the FLASH calculations are limited
            tmin_ht1 = max([multiRP.limitk(icomp=(i+1))['tmin'] \
                            for i in range(prop_hot_in['nc'])])
            #check limits
            tmin_ht2 = multiRP.limits(x_hot_in)['tmin']
        tmin_hot = max(tmin_ht1, tmin_ht2) + dt

        #check if multiprocessing is allowed
        if SetMultiprocessing().__repr__() == 'on':
            ##start with multiprocessing##
            if mRP == None:
                #create if not exist
                mRP = multiRP.multirefprop()#initiate multirefprop

            #create children
            de1 = mRP['process'](target=_dE1, args=(Q_cold_in, p_cold_in,
                                                    dp_cold, t_hot_in, dt,
                                                    tmax_cold, x_cold_in,
                                                    h_cold_in, prop_cold_in,
                                                    mRP))
            de2 = mRP['process'](target=_dE2, args=(q_cold_in, Q_cold_in,
                                                    Q_hot_in, p_cold_in,
                                                    p_hot_in, x_cold_in,
                                                    x_hot_in, h_cold_in,
                                                    h_hot_in, tmin_hot,
                                                    prop_cold_in, prop_hot_in,
                                                    dp_cold, dp_hot, mRP))
            de4 = mRP['process'](target=_dE4, args=(dp_hot, dt, h_hot_in,
                                                    p_hot_in, prop_hot_in,
                                                    Q_hot_in, t_cold_in,
                                                    tmin_hot, x_hot_in, mRP))
            de6 = mRP['process'](target=_dE6, args=(h_cold_in, h_hot_in,
                                                    p_cold_in, p_hot_in,
                                                    prop_cold_in, prop_hot_in,
                                                    Q_cold_in,q_hot_in,
                                                    Q_hot_in, tmax_cold,
                                                    x_cold_in, x_hot_in,
                                                    dp_cold, dp_hot, mRP))
            de7 = mRP['process'](target=_dE7, args=(q_cold_in, Q_cold_in,
                                                    Q_hot_in, p_cold_in,
                                                    p_hot_in, x_cold_in,
                                                    x_hot_in, h_cold_in,
                                                    h_hot_in, tmin_hot,
                                                    prop_cold_in, prop_hot_in,
                                                    dp_cold, dp_hot, mRP))
            de8 = mRP['process'](target=_dE8, args=(h_cold_in, h_hot_in,
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
            dE1 = mRP['result'][de1.name]
            dE2 = mRP['result'][de2.name]
            dE4 = mRP['result'][de4.name]
            dE6 = mRP['result'][de6.name]
            dE7 = mRP['result'][de7.name]
            dE8 = mRP['result'][de8.name]
            ##end with multiprocessing##

        #check if single processing is allowed
        elif SetMultiprocessing().__repr__() == 'off':
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
        if prop_in['t'] < prop_contra_in['t']:
            h_out = h_in + (dE / Q_cold_in)
            multiRP.resetup(prop_cold_in)
            prop = multiRP.flsh('ph', p_out, h_out, x_cold_in)

        #hot stream
        elif prop_in['t'] > prop_contra_in['t']:
            h_out = h_in - (dE / Q_hot_in)
            multiRP.resetup(prop_hot_in)
            prop = multiRP.flsh('ph', p_out, h_out, x_hot_in)

    #restore Q value
    prop['Q'] = Q_in

    #set error margin
    if prop_in['t'] < prop_contra_in['t']:
        if prop['Q'] == scipy.inf:
            E_cold = scipy.inf
        else:
            E_cold = abs(prop_in['h'] - prop['h']) * prop['Q']
        if prop_contra_in['Q'] == scipy.inf:
            E_hot = scipy.inf
        else:
            E_hot = abs(prop_contra_in['h'] -
                        prop_contra_out['h']) * prop_contra_in['Q']
    else:
        if prop['Q'] == scipy.inf:
            E_hot = scipy.inf
        else:
            E_hot = abs(prop_in['h'] - prop['h']) * prop['Q']
        if prop_contra_in['Q'] == scipy.inf:
            E_cold = scipy.inf
        else:
            E_cold = abs(prop_contra_in['h'] -
                         prop_contra_out['h']) * prop_contra_in['Q']
    #correction to prevent div/0 error
    if E_cold == 0 or E_hot == 0:
        errvalue = 1
    elif E_cold == scipy.inf or E_hot == scipy.inf:
        errvalue = 0
    elif E_cold >= E_hot:
        errvalue = 1 - abs(E_hot / E_cold)
    else:
        errvalue = 1 - abs(E_cold / E_hot)

    equipspec ={'type': 'EXCHANGER', 'name': name, 'dp': dp,
                'dp_contra': dp_contra, 'dt': dt, 'fluidin': prop_in,
                'fluidin_contra': prop_contra_in, 'fluidout': prop,
                'fluidout_contra': prop_contra_out, 'errvalue': errvalue,
                'E_hot':E_hot, 'E_cold':E_cold}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop

def _dE1(Q_cold_in, p_cold_in, dp_cold, t_hot_in, dt, tmax_cold, x_cold_in,
         h_cold_in, prop_cold_in, mRP):
    """define dE1 based on dt cold"""
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
        if t_cold_out >= crit['tcrit'] or p_cold_out >= crit['pcrit']:
            #flash calculation
            h_cold_out = multiRP.flsh('tp', t_cold_out,
                                      p_cold_out,
                                      x_cold_in)['h']
        else:
            #liquid condition
            if t_cold_out <= multiRP.satp(p_cold_out, x_cold_in, 1)['t']:
                D = multiRP.tprho(t_cold_out, p_cold_out, x_cold_in, 1)['D']
                h_cold_out = multiRP.enthal(t_cold_out, D, x_cold_in)['h']
            #vapor condition
            elif t_cold_out >= multiRP.satp(p_cold_out, x_cold_in, 2)['t']:
                D = multiRP.tprho(t_cold_out, p_cold_out, x_cold_in, 2)['D']
                h_cold_out = multiRP.enthal(t_cold_out, D, x_cold_in)['h']
            #two phase condition
            else:
                try:
                    #flash calculation
                    h_cold_out = multiRP.flsh('tp', t_cold_out,
                                              p_cold_out, x_cold_in)['h']
                except:
                    #asume liquid condition
                    ############################################################
                    # -1e-05 added due to round off errors at sat. bounderies. #
                    ############################################################
                    if t_cold_out - 1e-05 <= multiRP.satp(p_cold_out,
                                                          x_cold_in, 1)['t']:
                        D = multiRP.tprho(t_cold_out, p_cold_out,
                                          x_cold_in, 1)['D']
                        h_cold_out = multiRP.enthal(t_cold_out,
                                                    D, x_cold_in)['h']
                    else: raise EquipmentError()
        #heat exchanged
        dE1 = (h_cold_out - h_cold_in) * Q_cold_in
    #return results
    if SetMultiprocessing().__repr__() == 'on':
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = dE1
    elif SetMultiprocessing().__repr__() == 'off':
        return dE1

def _dE2(q_cold_in, Q_cold_in, Q_hot_in, p_cold_in, p_hot_in, x_cold_in,
         x_hot_in, h_cold_in, h_hot_in, tmin_hot, prop_cold_in, prop_hot_in,
         dp_cold, dp_hot, mRP):
    """define dE2 based on bubble point cold correction"""
    #settup cold fluid and initiate mRP
    multiRP.resetup(prop_cold_in, mRP=mRP)
    #check fluid not liquid
    phase = multiRP.getphase(prop_cold_in)
    if phase != 'liquid' and phase != 'saturated liquid':
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
        t_x_ing = prop['t']
        #check for cross temperature below min temp.
        if t_x_ing < tmin_hot:
            dE2 = scipy.inf
        else:
            #calculate enthalpy at crossing (cold)
            h = multiRP.therm(t_x_ing, prop['Dliq'], x_cold_in)['h']
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
                if t_x_ing >= crit['tcrit'] or p_hot_av >= crit['pcrit']:
                    #flash calculation
                    h = multiRP.flsh('tp', t_x_ing,
                                     p_hot_av, x_hot_in)['h']
                else:
                    #liquid condition
                    if t_x_ing <= multiRP.satp(p_hot_av, x_hot_in, 1)['t']:
                        D = multiRP.tprho(t_x_ing, p_hot_av, x_hot_in, 1)['D']
                        h = multiRP.enthal(t_x_ing, D, x_hot_in)['h']
                    #vapor condition
                    elif t_x_ing >= multiRP.satp(p_hot_av, x_hot_in, 2)['t']:
                        D = multiRP.tprho(t_x_ing, p_hot_av, x_hot_in, 2)['D']
                        h = multiRP.enthal(t_x_ing, D, x_hot_in)['h']
                    #two phase condition
                    else:
                        #flash calculation
                        h = multiRP.flsh('tp', t_x_ing,
                                         p_hot_av, x_hot_in)['h']
                #dE of hot part
                dE2b = (h_hot_in - h) * Q_hot_in
                #check infiniti
                if dE2b < 0: dE2 = scipy.inf
                #sum dE hot and cold
                else: dE2 = dE2a + dE2b
    #return results
    if SetMultiprocessing().__repr__() == 'on':
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = dE2
    elif SetMultiprocessing().__repr__() == 'off':
        return dE2

def _dE4(dp_hot, dt, h_hot_in, p_hot_in, prop_hot_in, Q_hot_in, t_cold_in,
         tmin_hot, x_hot_in, mRP):
    """define dE4 based on dt hot"""
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
        if t_hot_out >= crit['tcrit'] or p_hot_out >= crit['pcrit']:
            #flash calculation
            h_hot_out = multiRP.flsh('tp', t_hot_out, p_hot_out, x_hot_in)['h']
        else:
            #liquid condition
            ############################################################
            # -1e-07 added due to round off errors at sat. bounderies. #
            ############################################################
            if t_hot_out - 1e-07 <= multiRP.satp(p_hot_out, x_hot_in, 1)['t']:
                D = multiRP.tprho(t_hot_out, p_hot_out, x_hot_in, 1)['D']
                h_hot_out = multiRP.enthal(t_hot_out, D, x_hot_in)['h']
            #vapor condition
            elif t_hot_out >= multiRP.satp(p_hot_out, x_hot_in, 2)['t']:
                D = multiRP.tprho(t_hot_out, p_hot_out, x_hot_in, 2)['D']
                h_hot_out = multiRP.enthal(t_hot_out, D, x_hot_in)['h']
            #two phase condition
            else:
                try:
                    #flash calculation
                    h_hot_out = multiRP.flsh('tp', t_hot_out, p_hot_out,
                                             x_hot_in)['h']
                except:
                    #asume liquid condition
                    ############################################################
                    # -1e-05 added due to round off errors at sat. bounderies. #
                    ############################################################
                    if t_hot_out - 1e-05 <= multiRP.satp(p_hot_out, x_hot_in,
                                                         1)['t']:
                        D = multiRP.tprho(t_hot_out, p_hot_out, x_hot_in,
                                          1)['D']
                        h_hot_out = multiRP.enthal(t_hot_out, D, x_hot_in)['h']
                    else: raise EquipmentError()
        #heat exchanged
        dE4 = (h_hot_in - h_hot_out) * Q_hot_in
    #return values
    if SetMultiprocessing().__repr__() == 'on':
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = dE4
    elif SetMultiprocessing().__repr__() == 'off':
        return dE4

def _dE6(h_cold_in, h_hot_in, p_cold_in, p_hot_in, prop_cold_in, prop_hot_in,
         Q_cold_in, q_hot_in, Q_hot_in, tmax_cold, x_cold_in, x_hot_in, dp_cold,
         dp_hot, mRP):
    """define dE6 based on dew point hot correction"""
    #setup hot fluid and initiate mRP
    multiRP.resetup(prop_hot_in, mRP=mRP)
    #check for q < 1
    phase = multiRP.getphase(prop_hot_in)
    if phase != "vapor" and phase != "saturated vapor" and phase != "gas":
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
        t_x_ing = prop['t']
        #check for temp > max allowed
        if t_x_ing > tmax_cold:
            dE6 = scipy.inf
        else:
            #enthalpy at crossing
            h = multiRP.therm(t_x_ing, prop['Dvap'], x_hot_in)['h']
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
                if t_x_ing >= crit['tcrit'] or p_cold_av >= crit['pcrit']:
                    #flash calculation
                    h = multiRP.flsh('tp', t_x_ing, p_cold_av, x_cold_in)['h']
                else:
                    #liquid condition
                    if t_x_ing <= multiRP.satp(p_cold_av, x_cold_in, 1)['t']:
                        D = multiRP.tprho(t_x_ing, p_cold_av, x_cold_in, 1)['D']
                        h = multiRP.enthal(t_x_ing, D, x_cold_in)['h']
                    #vapor condition
                    elif t_x_ing >= multiRP.satp(p_cold_av, x_cold_in, 2)['t']:
                        D = multiRP.tprho(t_x_ing, p_cold_av, x_cold_in, 2)['D']
                        h = multiRP.enthal(t_x_ing, D, x_cold_in)['h']
                    #two phase condition
                    else:
                        #flash calculation
                        h = multiRP.flsh('tp', t_x_ing, p_cold_av,
                                         x_cold_in)['h']
                #dE cold fluid
                dE6b = (h - h_cold_in) * Q_cold_in
                #check for dE < 0
                if dE6b < 0: dE6 = scipy.inf
                #calculate dE
                else: dE6 = dE6a + dE6b
    #return results
    if SetMultiprocessing().__repr__() == 'on':
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = dE6
    elif SetMultiprocessing().__repr__() == 'off':
        return dE6

def _dE7(q_cold_in, Q_cold_in, Q_hot_in, p_cold_in, p_hot_in, x_cold_in,
         x_hot_in, h_cold_in, h_hot_in, tmin_hot, prop_cold_in, prop_hot_in,
         dp_cold, dp_hot, mRP):
    """define dE7 based on critical point cold correction"""
    #settup cold fluid and initiate mRP
    multiRP.resetup(prop_cold_in, mRP=mRP)
    #check fluid above pcrit.
    phase = multiRP.getphase(prop_cold_in)
    if phase != "compressible liquid":
        dE7 = scipy.inf
    #check infiniti
    elif Q_cold_in == scipy.inf or Q_hot_in == scipy.inf:
        dE7 = scipy.inf
    #calculate dE
    else:
        #calculate critical cond.
        prop = multiRP.critp(x_cold_in)
        #temp at potential crossing of two fluids
        t_x_ing = prop['tcrit']
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
            h = multiRP.flsh('tp', t_x_ing, p_cold_av, x_cold_in)['h']
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
                if t_x_ing >= crit['tcrit'] or p_hot_av >= crit['pcrit']:
                    #flash calculation
                    h = multiRP.flsh('tp', t_x_ing, p_hot_av, x_hot_in)['h']
                else:
                    #liquid condition
                    if t_x_ing <= multiRP.satp(p_hot_av, x_hot_in, 1)['t']:
                        D = multiRP.tprho(t_x_ing, p_hot_av, x_hot_in, 1)['D']
                        h = multiRP.enthal(t_x_ing, D, x_hot_in)['h']
                    #vapor condition
                    elif t_x_ing >= multiRP.satp(p_hot_av, x_hot_in, 2)['t']:
                        D = multiRP.tprho(t_x_ing, p_hot_av, x_hot_in, 2)['D']
                        h = multiRP.enthal(t_x_ing, D, x_hot_in)['h']
                    #two phase condition
                    else:
                        #flash calculation
                        h = multiRP.flsh('tp', t_x_ing, p_hot_av, x_hot_in)['h']
                #dE of hot part
                dE7b = (h_hot_in - h) * Q_hot_in
                #check infiniti
                if dE7b < 0: dE7 = scipy.inf
                #sum dE hot and cold
                else: dE7 = dE7a + dE7b
    #return results
    if SetMultiprocessing().__repr__() == 'on':
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = dE7
    elif SetMultiprocessing().__repr__() == 'off':
        return dE7

def _dE8(h_cold_in, h_hot_in, p_cold_in, p_hot_in, prop_cold_in, prop_hot_in,
         Q_cold_in, q_hot_in, Q_hot_in, tmax_cold, x_cold_in, x_hot_in, dp_cold,
         dp_hot, mRP):
    """define dE8 based on critical point hot correction"""
    #setup hot fluid and initiate mRP
    multiRP.resetup(prop_hot_in, mRP=mRP)
    #check for q < 1
    phase = multiRP.getphase(prop_hot_in)
    if phase != "Supercritical fluid":
        dE8 = scipy.inf
    #check for infinity
    elif Q_cold_in == scipy.inf or Q_hot_in == scipy.inf:
        dE8 = scipy.inf
    else:
        #calculate critical point
        prop = multiRP.critp(x_hot_in)
        #potential crossing temp
        t_x_ing = prop['tcrit']
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
            h = multiRPflsh('tp', t_x_ing, p_hot_av, x_hot_in)['h']
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
                if t_x_ing >= crit['tcrit'] or p_cold_av >= crit['pcrit']:
                    #flash calculation
                    h = multiRP.flsh('tp', t_x_ing, p_cold_av, x_cold_in)['h']
                else:
                    #liquid condition
                    if t_x_ing <= multiRP.satp(p_cold_av, x_cold_in, 1)['t']:
                        D = multiRP.tprho(t_x_ing, p_cold_av, x_cold_in, 1)['D']
                        h = multiRP.enthal(t_x_ing, D, x_cold_in)['h']
                    #vapor condition
                    elif t_x_ing >= multiRP.satp(p_cold_av, x_cold_in, 2)['t']:
                        D = multiRP.tprho(t_x_ing, p_cold_av, x_cold_in, 2)['D']
                        h = multiRP.enthal(t_x_ing, D, x_cold_in)['h']
                    #two phase condition
                    else:
                        #flash calculation
                        h = multiRP.flsh('tp', t_x_ing, p_cold_av,
                                         x_cold_in)['h']
                #dE cold fluid
                dE8b = (h - h_cold_in) * Q_cold_in
                #check for dE < 0
                if dE8b < 0: dE8 = scipy.inf
                #calculate dE
                else: dE8 = dE8a + dE8b
    #return results
    if SetMultiprocessing().__repr__() == 'on':
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = dE8
    elif SetMultiprocessing().__repr__() == 'off':
        return dE8


def turbine_separator(prop_in, p_out, dp, eff=1.0, e_eff=1.0, wet=0.0,
                      name=False, mRP=None):
    """Calculate fluid properties at turbine_separator outlet
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
    _inputerrorcheck(locals(), 'turbine_separator')
    multiRP.resetup(prop_in, mRP=mRP)

    #separator
    prop_liq, prop_vap = separator(prop_in, dp, mRP=mRP)

    ##start with multiprocessing##
    if SetMultiprocessing().__repr__() == 'on':
        if mRP == None:
            #create if not exist
            mRP = multiRP.multirefprop()#initiate multirefprop

        #create children
        tursep_reg = mRP['process'](target=reg_valve, args=(prop_liq, p_out),
                                    kwargs={'mRP':mRP})
        tursep_tur = mRP['process'](target=turbine, args=(prop_vap, p_out),
                                    kwargs={'mRP':mRP, 'eff':eff, 'e_eff':e_eff,
                                            'wet':wet})
        tursep_list = [tursep_reg, tursep_tur]

        #start and join children
        multiRP.run_mRP(tursep_list)

        #assign results
        prop_reg = mRP['result'][tursep_reg.name]
        prop_tur = mRP['result'][tursep_tur.name]
        ##end with multiprocessing##

    #or if single core
    elif SetMultiprocessing().__repr__() == 'off':
        prop_reg = reg_valve(prop_liq, p_out, mRP=None)
        prop_tur = turbine(prop_vap, p_out, mRP=None, eff=eff, e_eff=e_eff,
                           wet=wet)

    #set increased back pressure
    if 'backpress' in prop_tur:
        incr_back_press = prop_tur['backpress']
        prop_tur.pop('backpress')

    #flowmerge
    if prop_tur['Q'] == 0:
        prop = prop_reg
    elif prop_reg['Q'] == 0:
        prop = prop_tur
    else:
        prop = flowmerge(prop_reg, prop_tur)

    #power generated
    pwr = ((prop_vap['h'] - prop_tur['h']) * prop_vap['Q']) * e_eff

    #create equipment spec
    equipspec ={'type': 'TURBINESEP', 'name': name, 'eff': eff, 'p_out': p_out,
                'fluidin': prop_in, 'fluidout': prop, 'pwr': pwr, 'dp': dp,
                'e_eff':e_eff, 'wet':wet, 'incr_back_press':incr_back_press}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop


def turbine_bleed_sep(prop_in, p_out, P_outbleed, dp, Q_ratio, Q_ratiobleed,
                      eff=1.0, e_eff=1.0, wet=0.0, name=False, mRP=None):
    """Calculate fluid properties at turbine_separator outlet
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
    _inputerrorcheck(locals(), 'turbine_bleed_sep')
    multiRP.resetup(prop_in, mRP=mRP)

    #separator
    prop_liq, prop_vap = separator(prop_in, dp)

    #check if fluid in only liquid
    if prop_vap['Q'] == 0:
        ##start with multiprocessing##
        if SetMultiprocessing().__repr__() == 'on':
            if mRP == None:
                #create if not exist
                mRP = multiRP.multirefprop()#initiate multirefprop

            #create children
            tursep_0vap = mRP['process'](target=reg_valve,
                                         args=(prop_liq, p_out),
                                         kwargs={'mRP':mRP})
            tursep_0vapbleed = mRP['process'](target=reg_valve,
                                              args=(prop_liq, p_outbleed),
                                              kwargs={'mRP':mRP})
            tursep_list = [tursep_0vap, tursep_0vapbleed]

            #start and join children
            multiRP.run_mRP(tursep_list)

            #assign results
            prop = mRP['result'][tursep_0vap.name]
            prop_bleed = mRP['result'][tursep_0vapbleed.name]
            ##end with multiprocessing##

        #single core
        elif SetMultiprocessing().__repr__() == 'off':
            prop = reg_valve(prop_liq, p_out, mRP=None)
            prop_bleed = reg_valve(prop_liq, p_outbleed, mRP=None)

        #correct flow value
        prop['Q'] *= 1 - Q_ratiobleed
        prop_bleed['Q'] *= Q_ratiobleed

    #check if fluid contains vapor
    else:
        #copy prop_vap to bleed properties and modify Q value
        prop_bleed = copy.copy(prop_vap)
        prop_vap['Q'] *= 1 - Q_ratio
        prop_bleed['Q'] *= Q_ratio
        #copy prop_liq to bleed properties and modify Q value
        prop_liq_bleed = copy.copy(prop_liq)
        prop_liq['Q'] *= 1 - Q_ratiobleed
        prop_liq_bleed['Q'] *= Q_ratiobleed

        ##start with multiprocessing##
        if SetMultiprocessing().__repr__() == 'on':
            if mRP == None:
                #create if not exist
                mRP = multiRP.multirefprop()#initiate multirefprop

            #create children
            tursep_reg = mRP['process'](target=reg_valve,
                                        args=(prop_liq, p_out),
                                        kwargs={'mRP':mRP})
            tursep_reg_bleed = mRP['process'](target=reg_valve,
                                              args=(prop_liq_bleed, p_outbleed),
                                              kwargs={'mRP':mRP})
            tursep_tur = mRP['process'](target=turbine,
                                        args=(prop_vap, p_out),
                                        kwargs={'mRP':mRP, 'eff':eff,
                                                'e_eff':e_eff, 'wet':wet})
            tursep_bleed = mRP['process'](target=turbine,
                                          args=(prop_bleed, p_outbleed),
                                          kwargs={'mRP':mRP, 'eff':eff,
                                                  'e_eff':e_eff, 'wet':wet})
            tursep_list = [tursep_reg, tursep_reg_bleed, tursep_tur,
                           tursep_bleed]

            #start and join children
            multiRP.run_mRP(tursep_list)

            #assign results
            prop_reg = mRP['result'][tursep_reg.name]
            prop_reg_bleed = mRP['result'][tursep_reg_bleed.name]
            prop_tur = mRP['result'][tursep_tur.name]
            prop_bleed = mRP['result'][tursep_bleed.name]
            ##end with multiprocessing##

        #single core
        elif SetMultiprocessing().__repr__() == 'off':
            prop_reg = reg_valve(prop_liq, p_out, mRP=None)
            prop_reg_bleed = reg_valve(prop_liq_bleed, p_outbleed, mRP=None)
            prop_tur = turbine(prop_vap, p_out, mRP=None, eff=eff, e_eff=e_eff,
                               wet=wet)
            prop_bleed = turbine(prop_bleed, p_outbleed, mRP=None, eff=eff,
                                 e_eff=e_eff, wet=wet)

        #set increased back pressure
        incr_back_press = prop_tur['backpress']
        incr_back_press_bleed = prop_bleed['backpress']

        #flowmerge
        ##start with multiprocessing##
        if SetMultiprocessing().__repr__() == 'on':
            if mRP == None:
                #create if not exist
                mRP = multiRP.multirefprop()#initiate multirefprop

            #create children
            tursep_tur = mRP['process'](target=flowmerge,
                                        args=(prop_reg, prop_tur),
                                        kwargs={'mRP':mRP})
            tursep_bleed = mRP['process'](target=flowmerge,
                                          args=(prop_reg_bleed, prop_bleed),
                                          kwargs={'mRP':mRP})
            tursep_list = [tursep_tur, tursep_bleed]

            #start and join children
            multiRP.run_mRP(tursep_list)

            #assign results
            prop = mRP['result'][tursep_tur.name]
            prop_bleed = mRP['result'][tursep_bleed.name]
            ##end multiprocessing##

        #single core
        elif SetMultiprocessing().__repr__() == 'off':
            prop = flowmerge(prop_reg, prop_tur, mRP=None)
            prop_bleed = flowmerge(prop_reg_bleed, prop_bleed, mRP=None)

    #power generated
    pwr = (((prop_vap['h'] - prop_tur['h']) * prop_vap['Q']) * e_eff +
          ((prop_vap['h'] - prop_bleed['h']) * prop_bleed['Q']) * e_eff)

    #create equipment spec
    equipspec ={'type': 'TURBINEBLEEDSEP', 'name': name, 'eff': eff,
                'p_out': p_out, 'fluidin': prop_in, 'fluidout': prop,
                'pwr': pwr, 'dp': dp, 'e_eff':e_eff, 'p_outbleed':p_outbleed,
                'Q_ratio':Q_ratio, 'fluidbleedout':prop_bleed,
                'Q_ratiobleed':Q_ratiobleed, 'wet':wet,
                'incr_back_press':incr_back_press,
                'incr_back_press_bleed':incr_back_press_bleed}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = (prop,
                                                                    prop_bleed)

    #return values
    if name != False:
        return prop, prop_bleed, equipspec
    else: return prop, prop_bleed


def reg_valve(prop_in, p_out, name=False, mRP=None):
    """Calculate fluid properties downstream press. reg. valve with constant
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
    _inputerrorcheck(locals(), 'reg_valve')
    multiRP.resetup(prop_in, mRP=mRP)
    h = prop_in['h']
    x = prop_in['x']

    #calculate properties downstream valve
    prop = multiRP.flsh('ph', p_out, h, x)

    #add Q
    if 'Q' in prop_in:
        prop['Q'] = prop_in['Q']

    equipspec ={'type': 'REGVALVE', 'name': name, 'p_out': p_out,
                'fluidin': prop_in, 'fluidout': prop}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop


def separator(prop_in, dp, name=False, mRP=None):
    """Calculate fluid properties at separator outlets (liquid and vapor).
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
    _inputerrorcheck(locals(), 'separator')
    multiRP.resetup(prop_in, mRP=mRP)
    p = prop_in['p'] - dp
    x = prop_in['x']
    h = prop_in['h']
    Q = prop_in['Q']

    #calculate properties with corrected p and constant h
    #enhance speed
    if dp > 0:
        #inhibit error reporting
        if str(multiRP.SetErrorDebug()) == 'on':
            sed = multiRP.SetErrorDebug.on
            multiRP.SetErrorDebug.off()
        elif str(multiRP.SetErrorDebug()) == 'off':
            sed = multiRP.SetErrorDebug.off
        try:
            prop = multiRP.ph2ph(p, h, x)
            prop.update(sed())
        except multiRP.RefpropError:
            sed() #re-activate error reporting
            prop = multiRP.flsh('ph', p, h, x)
    else: prop = copy.copy(prop_in)

    prop_liq = copy.copy(prop)
    prop_vap = copy.copy(prop)
    prop_liq.update({'D': prop['Dliq'], 'x': prop['xliq']})
    prop_vap.update({'D': prop['Dvap'], 'x': prop['xvap']})

    #vapor quality to be between 0 and 1
    phase = multiRP.getphase(prop)
    if phase == "liquid" or phase == "compressible liquid" \
    or phase == "saturated liquid":
        q = 0
    elif phase == "saturated vapor" or phase == "gas" \
    or phase == "Supercritical fluid" or phase == "vapor":
        q = 1
    elif phase == "2 phase":
        q = prop['q']

    ##start with multiprocessing##
    if SetMultiprocessing().__repr__() == 'on':
        if mRP == None:
            #create mRP if not exist
            mRP = multiRP.multirefprop()#initiate multirefprop

        #create children
        sep_liq = mRP['process'](target=_separator_liq, args=(prop_liq, p, mRP))
        sep_vap = mRP['process'](target=_separator_vap, args=(prop_vap, p, mRP))
        sep_list = [sep_liq, sep_vap]

        #start and join children
        multiRP.run_mRP(sep_list)

        #assign results
        prop_liq = mRP['result'][sep_liq.name]
        prop_vap = mRP['result'][sep_vap.name]
        ##end with multiprocessing##

    #single core
    elif SetMultiprocessing().__repr__() == 'off':
        prop_liq = _separator_liq(prop_liq, p, None)
        prop_vap = _separator_vap(prop_vap, p, None)

    #add Q to input values for liquid and vapor
    prop_liq['Q'] = Q * (1 - q)
    prop_vap['Q'] = Q * q

    #create equipspec
    equipspec ={'type': 'SEPARATOR', 'name': name, 'dp': dp,
                'fluidin': prop_in, 'fluidliq': prop_liq, 'fluidvap': prop_vap}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = prop

    #return values
    if name != False:
        return prop_liq, prop_vap, equipspec
    else: return prop_liq, prop_vap

#child1 for separator, calculate liquid properties
def _separator_liq(prop_liq, p, mRP):
    'child calculation for def "separator"'
    prop = multiRP.therm(prop_liq['t'], prop_liq['D'], prop_liq['x'],
                         prop=prop_liq, mRP=mRP)
    prop['p'] = p
    prop['q'] = 0 #to saturated liquid state
    if SetMultiprocessing().__repr__() == 'on':
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = prop
    elif SetMultiprocessing().__repr__() == 'off':
        return prop

#child2 for separator, calculate vapor properties
def _separator_vap(prop_vap, p, mRP):
    'child calculation for def "separator"'
    prop = multiRP.therm(prop_vap['t'], prop_vap['D'], prop_vap['x'],
                         prop=prop_vap, mRP=mRP)
    prop['p'] = p
    prop['q'] = +1 #to saturated vapor state
    if SetMultiprocessing().__repr__() == 'on':
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = prop
    elif SetMultiprocessing().__repr__() == 'off':
        return prop


def compressor(prop_in, dp, eff=1.0, e_eff=1.0, name=False, mRP=None):
    """Calculate fluid properties at compressor discharge.
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
    _inputerrorcheck(locals(), 'compressor')
    multiRP.resetup(prop_in, mRP=mRP)
    p = prop_in['p'] + dp
    s = prop_in['s']
    x = prop_in['x']
    flshcalc = False
    if 'h' in prop_in: h = prop_in['h']
    else:
        #inhibit error reporting
        if str(multiRP.SetErrorDebug()) == 'on':
            sed = multiRP.SetErrorDebug.on
            multiRP.SetErrorDebug.off()
        elif str(multiRP.SetErrorDebug()) == 'off':
            sed = multiRP.SetErrorDebug.off
        try: #use flsh1 calculation for speed increase
            h = prop = multiRP.psvap(prop_in['p'], s, x)['h']
            prop.update(sed())
        except multiRP.RefpropError:
            sed() #re-activate error reporting
            flshcalc = True
            h = multiRP.flsh('ps', p, s, x)['h']
    Q = prop_in['Q']

    #inhibit error reporting
    if str(multiRP.SetErrorDebug()) == 'on':
        sed = multiRP.SetErrorDebug.on
        multiRP.SetErrorDebug.off()
    elif str(multiRP.SetErrorDebug()) == 'off':
        sed = multiRP.SetErrorDebug.off
    #calculate enthalpy at 100% efficiency
    try: #use flsh1 calculation for speed increase
        if flshcalc: raise multiRP.RefpropError
        prop = multiRP.psvap(p, s, x)
        prop.update(sed())
    except multiRP.RefpropError:
        sed() #re-activate error reporting
        flshcalc = True
        prop = multiRP.flsh('ps', p, s, x)

    #correct h for efficiency losses
    h_out = h + ((prop['h'] - h) / eff)

    #calculate properties with corrected h and discharge p
    #inhibit error reporting
    if str(multiRP.SetErrorDebug()) == 'on':
        sed = multiRP.SetErrorDebug.on
        multiRP.SetErrorDebug.off()
    elif str(multiRP.SetErrorDebug()) == 'off':
        sed = multiRP.SetErrorDebug.off
    try: #use flsh1 calculation for speed increase
        if flshcalc: raise multiRP.RefpropError
        prop = multiRP.phvap(p, h_out, x)
        prop.update(sed())
    except multiRP.RefpropError:
        sed() #re-activate error reporting
        prop = multiRP.flsh('ph', p, h_out, x)

    #restore Q value
    prop['Q'] = Q

    #power consumption
    pwr = ((h_out - h) * Q) / e_eff

    equipspec ={'type': 'COMPRESSOR', 'name': name, 'dp': dp, 'eff': eff,
                'fluidin': prop_in, 'fluidout': prop, 'pwr': -pwr,
                'e_eff':e_eff}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop



def subcool(prop_in, dt, name=False, mRP=None):
    """Calculate fluid properties based on subcooled temperature

    input:
        prop--fluid properties library contain keys:
            't'--subcool temperature.
            'x'--fluid composition
        dt--temperature difference between bubble point and subcool
    output:
        prop--fluid out properties
        equipspec--equipment details"""
    _inputerrorcheck(locals(), 'subcool')
    multiRP.resetup(prop_in, mRP=mRP)

    #calculate saturated properties at elevated temperature
    p = multiRP.satt(prop_in['t'] + dt, prop_in['x'], 1)['p']

    #calculate sub cooled properties at set temp and calc. press
    try:
        D = multiRP.tprho(prop_in['t'], p, prop_in['x'], 1)['D']
        prop = multiRP.therm(prop_in['t'], D, prop_in['x'])
        prop['q'] = -1
        prop['p'] = p
    except:
        prop = multiRP.flsh('tp', prop_in['t'], p, prop_in['x'])

    #add Q
    if 'Q' in prop_in:
        prop['Q'] = prop_in['Q']

    #create equipment spec
    equipspec ={'type': 'SUBCOOL', 'name': name, 'dt': dt, 'fluidin': prop_in,
                'fluidout': prop}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop


def turbine(prop_in, p_out, eff=1.0, e_eff=1.0, wet=0.0, name=False, mRP=None):
    """Calculate fluid properties at turbine discharge.
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
    _inputerrorcheck(locals(), 'turbine')
    multiRP.resetup(prop_in, mRP=mRP)
    p = p_out
    s = prop_in['s']
    x = prop_in['x']
    h = prop_in['h']
    Q = prop_in['Q']

    def _turbine():
        #set error
        flsh_2ph = False
        flsh_vap = False
        #calculate enthalpy at 100% efficiency
        #inhibit error reporting
        if str(multiRP.SetErrorDebug()) == 'on':
            sed = multiRP.SetErrorDebug.on
            multiRP.SetErrorDebug.off()
        elif str(multiRP.SetErrorDebug()) == 'off':
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
                prop = multiRP.flsh('ps', p, s, x)
        #correct h for efficiency losses
        h_out = h - ((h - prop['h']) * eff)

        #calculate properties with corrected h and discharge p
        #inhibit error reporting
        if str(multiRP.SetErrorDebug()) == 'on':
            sed = multiRP.SetErrorDebug.on
            multiRP.SetErrorDebug.off()
        elif str(multiRP.SetErrorDebug()) == 'off':
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
                prop = multiRP.flsh('ph', p, h_out, x)
        #restore input values
        prop['Q'] = Q
        #power generated
        pwr = ((h - h_out) * prop_in['Q']) * e_eff
        #return
        return (prop, pwr)

    #run turbine
    prop, pwr = _turbine()

    #check wet gas within limits
    if prop['q'] >= wet:
        incr_back_press = 0
    else:
        #set initial press adjustment value
        p_adj = 1
        #set initial delta pressure
        press_delta = prop_in['p'] - p
        #create while loop to find optimum back pressure
        while not wet - 0.00001 < prop['q'] < wet + 0.00001:
            #adjust press_adjust with factor 2
            p_adj *= 2
            #store previous run
            q_old = prop['q']
            p_old = p
            #correct pressure delta value
            if prop['q'] <= wet:
                press_delta -= press_delta / p_adj
            elif prop['q'] >= wet:
                press_delta += press_delta / p_adj
            #new pressure value
            p = prop_in['p'] - press_delta
            #run turbine
            prop, pwr = _turbine()
            #compare new values against old to determin insoluble situation
            if q_old > prop['q'] and p_old < p:
                raise EquipmentinfeasibleError('error raised insoluble q value')
            if p_adj >= 2 ** 16:
                raise EquipmentinfeasibleError('error raised unsolvable q value')
        #additional press
        incr_back_press = p - p_out
        #depress to reach p_out
        prop = reg_valve(prop, p_out, mRP=mRP)

    #create equipment spec
    equipspec ={'type': 'TURBINE', 'name': name, 'eff': eff, 'p_out': p_out,
                'fluidin': prop_in, 'fluidout': prop, 'pwr': pwr,
                'e_eff':e_eff, 'wet':wet, 'incr_back_press':incr_back_press}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        prop['backpress'] = incr_back_press
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = prop
        prop.pop('backpress')

    #return values
    if name != False:
        return prop, equipspec
    else:
        #add back pressure for turbine separator
        prop['backpress'] = incr_back_press
        return prop


def flowmerge(prop_1, prop_2=None, name=False, mRP=None):
    """Calculate fluid properties of 2 flow merge.
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
    _inputerrorcheck(locals(), 'flowmerge')

    #return prop_1 value if prop_2 had been deleted
    if not 'prop_2' in locals():
        prop = copy.copy(prop_1)

        #create equipment spec
        equipspec ={'type': 'FLOWMERGE', 'name': name, 'fluid1': prop_1,
                    'fluid2': None, 'fluidout': prop}

    else:
        #confirm both prop input are same refprop setup
        if multiRP.setup_details(prop_1) != multiRP.setup_details(prop_2):
            raise EquipmentinputError('''prop_1 and prop_2 input are not
                                      matching, ensure input are of same
                                      setup''')

        multiRP.resetup(prop_1, mRP=mRP)
        p1 = prop_1['p']
        x1 = prop_1['x']
        h1 = decimal.Decimal(prop_1['h'])
        Q1 = decimal.Decimal(prop_1['Q'])
        p2 = prop_2['p']
        x2 = prop_2['x']
        h2 = decimal.Decimal(prop_2['h'])
        Q2 = decimal.Decimal(prop_2['Q'])
        Q = Q1 + Q2

        #calculate merged x value.
        x = [(x1[each] * Q1 + x2[each] * Q2) / Q for each in range(len(x1))]
        x = multiRP.normalize(x)['x']

        #calculate merged h value.
        h = float((h1 * Q1 + h2 * Q2) / Q)

        #calculate fluid properties
        prop = multiRP.flsh('ph', p1, h, x)

        #calculate Q
        prop['Q'] = float(Q)

        #create equipment spec
        equipspec ={'type': 'FLOWMERGE', 'name': name, 'fluid1': prop_1,
                    'fluid2': prop_2, 'fluidout': prop}

    #return value if multiprocessing and if child process
    if mRP != None and mp.current_process()._parent_pid != None:
        multiRP.result[multiRP.process().name.rsplit(':', 1)[0]] = prop

    #return values
    if name != False:
        return prop, equipspec
    else: return prop


def pwr_gen(equipspecs):
    """Returns the total power generated - the total power consumed

    input:
        equipspecs--equipspec output from equipment in list format"""
    pwr = 0
    for key in equipspecs:
        if 'pwr' in equipspecs[key].keys():
            pwr += equipspecs[key]['pwr']
    return pwr


def exch_err(equipspecs, criteria=None):
    """Returns error margin of exchangers:
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
        if equipspecs[key]['type'] == 'EXCHANGER':
            err.append(equipspecs[key]['errvalue'])
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
        raise EquipmentinputError("input value criteria to be between 0 and " +
                                    "1 instead of " + str(criteria))


def cycletotal(equipspecs, prescript='', postscript='', mRP=None):
    """returns a full overview of total power consumption and all
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
        if equipspecs[key]['type'] == 'TURBINE':
            turbines.append(equipspecs[key])
        elif equipspecs[key]['type'] == 'TURBINESEP':
            turbineseps.append(equipspecs[key])
        elif equipspecs[key]['type'] == 'TURBINEBLEEDSEP':
            turbinebleedseps.append(equipspecs[key])
        elif equipspecs[key]['type'] == 'PUMP':
            pumps.append(equipspecs[key])
        elif equipspecs[key]['type'] == 'COMPRESSOR':
            compressors.append(equipspecs[key])
        elif equipspecs[key]['type'] == 'EXCHANGER':
            exchangers.append(equipspecs[key])
        elif equipspecs[key]['type'] == 'SUBCOOL':
            subcools.append(equipspecs[key])
        elif equipspecs[key]['type'] == 'SEPARATOR':
            separators.append(equipspecs[key])
        elif equipspecs[key]['type'] == 'REGVALVE':
            regvalves.append(equipspecs[key])
        elif equipspecs[key]['type'] == 'FLOWMERGE':
            flowmerges.append(equipspecs[key])
    #equipment string
    equipstr = ''

    #power string
    powergen = 0
    powercon = 0
    pwrstr = '*' * 80
    pwrstr += '{:}{:^80}{:}'.format('\n', 'POWER GENERATION', '\n')

    #turbine
    equipdisplay = True
    for turbine in turbines:
        #power string cont.
        #ensure Turbine is displayed once and only if turbine type in equipspec
        if equipdisplay:
            pwrstr += 'TURBINE\n'
            equipdisplay = False
        pwrstr += _prop(turbine['name'], turbine['pwr'], '  Watt')
        powergen += turbine['pwr']
        #equipment string cont.
        equipstr += '*' * 80 + '\n'
        equipstr += '{:^80}'.format('TURBINE') + '\n'
        equipstr+= '{:^80}'.format(turbine['name']) + '\n\n'
        equipstr += _prop('Shaft power output: ', turbine['pwr'],
                          '  Watt')
        equipstr += _prop('Turbine eff.: ', turbine['eff'] * 100, '%')
        equipstr += _prop('Alter. eff.: ', turbine['e_eff'] * 100, '%')
        equipstr += _prop('Min. vapor outl.: ', turbine['wet'] * 100,
                          ' mol(v)/mol')
        equipstr += _prop('add. back press.: ', turbine['incr_back_press'],
                          ' kPa')
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID IN', '\n')
        equipstr += _fluidprop(turbine['fluidin'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID OUT', '\n')
        equipstr += _fluidprop(turbine['fluidout'])

    #turbine separator
    equipdisplay = True
    for turbinesep in turbineseps:
        #power string cont.
        #ensure Turbine is displayed once and only if turbine type in equipspec
        if equipdisplay:
            pwrstr += 'TURBINE SEPARATOR\n'
            equipdisplay = False
        pwrstr += _prop(turbinesep['name'], turbinesep['pwr'], '  Watt')
        powergen += turbinesep['pwr']
        #equipment string cont.
        equipstr += '*' * 80 + '\n'
        equipstr += '{:^80}'.format('TURBINE SEPARATOR') + '\n'
        equipstr+= '{:^80}'.format(turbinesep['name']) + '\n\n'
        equipstr += _prop('Shaft power output: ', turbinesep['pwr'],
                          '  Watt')
        equipstr += _prop('Turbine eff.: ', turbinesep['eff'] * 100, '%')
        equipstr += _prop('Alter. eff.: ', turbinesep['eff'] * 100, '%')
        equipstr += _prop('separator dp: ', turbinesep['dp'], ' kPa')
        equipstr += _prop('Min. vapor outl.: ', turbinesep['wet'] * 100,
                          ' mol(v)/mol')
        equipstr += _prop('add. back press.: ', turbinesep['incr_back_press'],
                          ' kPa')
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID IN', '\n')
        equipstr += _fluidprop(turbinesep['fluidin'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID OUT', '\n')
        equipstr += _fluidprop(turbinesep['fluidout'])

    #power string cont.
    pwrstr += _prop('Subtotal: ', powergen, '  Watt')
    pwrstr += '{:}{:^80}{:}'.format('\n', 'POWER CONSUMPTION', '\n')

    #turbine bleed separator
    equipdisplay = True
    for turbinebleedsep in turbinebleedseps:
        #power string cont.
        #ensure Turbine is displayed once and only if turbine type in equipspec
        if equipdisplay:
            pwrstr += 'TURBINE BLEED SEPARATOR\n'
            equipdisplay = False
        pwrstr += _prop(turbinebleedsep['name'], turbinebleedsep['pwr'],
                        '  Watt')
        powergen += turbinebleedsep['pwr']
        #equipment string cont.
        equipstr += '*' * 80 + '\n'
        equipstr += '{:^80}'.format('TURBINE BLEED SEPARATOR') + '\n'
        equipstr+= '{:^80}'.format(turbinebleedsep['name']) + '\n\n'
        equipstr += _prop('Shaft power output: ', turbinebleedsep['pwr'],
                          '  Watt')
        equipstr += _prop('Turbine eff.: ', turbinebleedsep['eff'] * 100, '%')
        equipstr += _prop('Alter. eff.: ', turbinebleedsep['eff'] * 100, '%')
        equipstr += _prop('separator dp: ', turbinebleedsep['dp'], ' kPa')
        equipstr += _prop('flow ratio: ',
                          turbinebleedsep['Q_ratio'] * 100, '%')
        equipstr += _prop('cond. flow ratio: ',
                          turbinebleedsep['Q_ratiobleed'] * 100, '%')
        equipstr += _prop('Min. vapor outl.: ', turbinebleedsep['wet'] * 100,
                          ' mol(v)/mol')
        equipstr += _prop('add. back press.: ',
                          turbinebleedsep['incr_back_press'], ' kPa')
        equipstr += _prop('add. bk pr. bleed: ',
                          turbinebleedsep['incr_back_press'], ' kPa')
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID IN', '\n')
        equipstr += _fluidprop(turbinebleedsep['fluidin'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID BLEED OUT', '\n')
        equipstr += _fluidprop(turbinebleedsep['fluidbleedout'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID OUT', '\n')
        equipstr += _fluidprop(turbinebleedsep['fluidout'])

    #power string cont.
    pwrstr += _prop('Subtotal: ', powergen, '  Watt')
    pwrstr += '{:}{:^80}{:}'.format('\n', 'POWER CONSUMPTION', '\n')

    #pump
    equipdisplay = True
    for pump in pumps:
        #power string cont.
        #ensure pump is displayed once and only if pump type in equipspec
        if equipdisplay:
            pwrstr += 'PUMP\n'
            equipdisplay = False
        pwrstr += _prop(pump['name'], -pump['pwr'], '  Watt')
        powercon += -pump['pwr']
        #equipment string cont.
        equipstr += '*' * 80 + '\n'
        equipstr += '{:^80}'.format('PUMP') + '\n'
        equipstr += '{:^80}'.format(pump['name']) + '\n\n'
        equipstr += _prop('Power consumed: ', -pump['pwr'], '  Watt')
        equipstr += _prop('Pump eff.: ', pump['eff'] * 100, '%')
        equipstr += _prop('motor eff.: ', pump['e_eff'] * 100, '%')
        equipstr += _prop('pump head: ', pump['dp'], ' kPa', pump['dp'] / 100,
                          ' bar(a)')
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID IN', '\n')
        equipstr += _fluidprop(pump['fluidin'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID OUT', '\n')
        equipstr += _fluidprop(pump['fluidout'])

    #compressor
    equipdisplay = True
    for compressor in compressors:
        #power string cont.
        #ensure pump is displayed once and only if pump type in equipspec
        if equipdisplay:
            pwrstr += 'COMPRESSOR\n'
            equipdisplay = False
        pwrstr += _prop(compressor['name'], -compressor['pwr'], '  Watt')
        powercon += -compressor['pwr']
        #equipment string cont.
        equipstr += '*' * 80 + '\n'
        equipstr += '{:^80}'.format('COMPRESSOR') + '\n'
        equipstr += '{:^80}'.format(compressor['name']) + '\n\n'
        equipstr += _prop('Power consumed: ', -compressor['pwr'], '  Watt')
        equipstr += _prop('Compr. eff.: ', compressor['eff'] * 100, '%')
        equipstr += _prop('Motor eff.: ', compressor['e_eff'] * 100, '%')
        equipstr += _prop('compressor head: ', compressor['dp'], ' kPa')
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID IN', '\n')
        equipstr += _fluidprop(compressor['fluidin'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID OUT', '\n')
        equipstr += _fluidprop(compressor['fluidout'])

    #power string cont.
    pwrstr += _prop('Subtotal: ', powercon, '  Watt')
    pwrnett = powergen - powercon
    if pwrnett == scipy.inf:
        pwrnett = -scipy.inf
    pwrstr += '{:}{:^80}{:}'.format('\n', 'NETT POWER GENERATION', '\n')
    pwrstr += _prop('Nett power gen.: ', pwrnett, '  Watt')

    #exchanger
    for exchanger in exchangers:
        #equipment string cont.
        equipstr += '*' * 80 + '\n'
        equipstr += '{:^80}'.format('EXCHANGER') + '\n'
        equipstr += '{:^80}'.format(exchanger['name']) + '\n\n'
        equipstr += _prop('Cooling duty (cold): ', exchanger['E_cold'],
                          '  Watt')
        equipstr +=_prop('Cooling duty (hot): ', exchanger['E_hot'],
                         '  Watt')
        equipstr += _prop('Error value: ', exchanger['errvalue'] * 100,
                          '%')
        equipstr += _prop('Delta temp.: ', exchanger['dt'], '  Kelvin')
        #if fluidin = cold
        if exchanger['fluidin']['t'] < exchanger['fluidin_contra']['t']:
            equipstr += _prop('Delta press. cold: ', exchanger['dp'], '  kPa')
            equipstr += _prop('Delta press. hot: ', exchanger['dp_contra'],
                              '  kPa')
            equipstr += '{:}{:^80}{:}'.format('\n', 'COLD FLUID IN', '\n')
            equipstr += _fluidprop(exchanger['fluidin'])
            equipstr += '{:}{:^80}{:}'.format('\n', 'COLD FLUID OUT', '\n')
            equipstr += _fluidprop(exchanger['fluidout'])
            equipstr += '{:}{:^80}{:}'.format('\n', 'HOT FLUID IN', '\n')
            equipstr += _fluidprop(exchanger['fluidin_contra'])
            equipstr += '{:}{:^80}{:}'.format('\n', 'HOT FLUID OUT', '\n')
            equipstr += _fluidprop(exchanger['fluidout_contra'])
        #if fluidin = hot
        else:
            equipstr += _prop('Delta press. cold: ', exchanger['dp_contra'],
                              '  kPa')
            equipstr += _prop('Delta press. hot: ', exchanger['dp'], '  kPa')
            equipstr += '{:}{:^80}{:}'.format('\n', 'COLD FLUID IN', '\n')
            equipstr += _fluidprop(exchanger['fluidin_contra'])
            equipstr += '{:}{:^80}{:}'.format('\n', 'COLD FLUID OUT', '\n')
            equipstr += _fluidprop(exchanger['fluidout_contra'])
            equipstr += '{:}{:^80}{:}'.format('\n', 'HOT FLUID IN', '\n')
            equipstr += _fluidprop(exchanger['fluidin'])
            equipstr += '{:}{:^80}{:}'.format('\n', 'HOT FLUID OUT', '\n')
            equipstr += _fluidprop(exchanger['fluidout'])
        #add graph
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLOW - TEMP.', '\n')
        equipstr += '{:^80}{:}'.format('GRAPHIC', '\n')
        try:
            graph = _exchgraph(exchanger['fluidin'],
                               exchanger['fluidout'],
                               exchanger['fluidin_contra'],
                               exchanger['fluidout_contra'], mRP=mRP)
        except:
            graph = 'error in generating graphics at def cycletotal\n'
        equipstr += graph

    #subcool
    for subcool in subcools:
        #equipment string cont.
        equipstr += '*' * 80 + '\n'
        equipstr += '{:^80}'.format('SUBCOOL') + '\n'
        equipstr += '{:^80}'.format(subcool['name']) + '\n\n'
        equipstr += _prop('Delta temp.: ', subcool['dt'], '  Kelvin')
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID IN', '\n')
        equipstr += _fluidprop(subcool['fluidin'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID OUT', '\n')
        equipstr += _fluidprop(subcool['fluidout'])

    #separator
    for separator in separators:
        #equipment strin cont.
        equipstr += '*' * 80 + '\n'
        equipstr += '{:^80}'.format('SEPARATOR') + '\n'
        equipstr += '{:^80}'.format(separator['name']) + '\n\n'
        equipstr += _prop('Delta press.: ', separator['dp'], '  kPa')
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID IN', '\n')
        equipstr += _fluidprop(separator['fluidin'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'VAPOR OUT', '\n')
        equipstr += _fluidprop(separator['fluidvap'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'LIQUID OUT', '\n')
        equipstr += _fluidprop(separator['fluidliq'])

    #regulating valve
    for regvalve in regvalves:
        #equipment strin cont.
        dp = regvalve['fluidin']['p'] - regvalve['p_out']
        equipstr += '*' * 80 + '\n'
        equipstr += '{:^80}'.format('REGULATING VALVE') + '\n'
        equipstr += '{:^80}'.format(regvalve['name']) + '\n\n'
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID IN', '\n')
        equipstr += _fluidprop(regvalve['fluidin'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID OUT', '\n')
        equipstr += _fluidprop(regvalve['fluidout'])

    #2-flow merge
    for flowmerge in flowmerges:
        #equipment strin cont.
        equipstr += '*' * 80 + '\n'
        equipstr += '{:^80}'.format('FLOWMERGE') + '\n'
        equipstr += '{:^80}'.format(flowmerge['name']) + '\n\n'
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID 1 IN', '\n')
        equipstr += _fluidprop(flowmerge['fluid1'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID 2 IN', '\n')
        equipstr += _fluidprop(flowmerge['fluid2'])
        equipstr += '{:}{:^80}{:}'.format('\n', 'FLUID OUT', '\n')
        equipstr += _fluidprop(flowmerge['fluidout'])

    #last line
    equipstr += '*' * 80 + '\n'

    #return total string
    return prescript + '\n' + pwrstr + equipstr + postscript


#for testing only
def _test(prop, mRP=None):

    equipdct = {}

    print('nett power before equipment generation')
    print(pwr_gen(equipdct))

    print('run subcool')
    prop, equipdct['subcool_01'] = subcool(prop, 4.5, 'subcool-01')

    print('run pump')
    prop03, equipdct['pump_01'] = pump(prop, 1000, 0.85, 0.94, name='pump-01')

    #exch_err(0.4)
    print('run exchanger simple')
    prop04, equipdct['exchanger_01'] = exchanger(prop03, 25, 4,
                                                 name='exchanger-01')

    print('run separator')
    prop04['t'] += 100
    prop01 = multiRP.flsh('tp', prop04['t'], prop04['p'], prop04['x'])
    prop01['Q'] = 1
    propliq, prop02, equipdct['separator_01'] = separator(prop01, 50,
                                                          'separator-01')

    print('run reg_valve')
    prop_reg_valve, equipdct['regvalve_01'] = reg_valve(propliq, 100,
                                                        'regvalve-01')

    print('run turbine')
    prop, equipdct['turbine_01'] = turbine(prop02, 100, 0.85, 0.96,
                                           wet=0.97, name='Turbine-01')

    print('run flowmerge')
    prop_fm, equipdct['flowmerge_01'] = flowmerge(prop, prop_reg_valve,
                                                  'flowmerge-01')

    print('run turbine_separator')
    prop, equipdct['turbine_separator_01'] = turbine_separator(prop01,
        100, 50, 0.85, wet=0.95, name='turbine_separator-01')

    print('run compressor & exchanger inletflow')
    prop05, equipdct['compressor_01']= compressor(prop02, 3000, 0.8,
                                                 0.94, name='compressor-01')

    print('run exchanger')
    prop06, equipdct['exchanger_01'] = exchanger(prop05, 25, 4, prop03, 25,
                                                 prop04, name='Exchanger-01')

    #exch_err(0.1)
    print('run exchanger otherflow outlet')
    prop07, equipdct['exchanger_01'] = exchanger(prop03, 25, 4, prop05, 25,
                                                 prop06, name='Exchanger-01')

    if abs((((prop05['h'] - prop06['h']) * prop05['Q']) -
             ((prop07['h'] - prop03['h']) * prop03['Q'])) /
            ((prop05['h'] - prop06['h']) * prop05['Q'])) < 0.00000001:
        print('exchanger tolerance test acceptable\n')
    else:
        print('exchanger tolerance test failed\n')

    print('exchanger err value')
    print(exch_err(equipdct))
    print(exch_err(equipdct, 0.02))


    print('nett power')
    print(pwr_gen(equipdct))

    return cycletotal(equipdct, prescript='start', postscript='end')

if __name__ == '__main__':
    #add module file path to python sys path
    import equipment as _filename
    _filename = (os.path.dirname(_filename.__file__))
    sys.path.append(_filename)

    #examples and test setup
    #setup fluid
    prop = multiRP.setup('def', 'water', 'ammonia')
    prop['x'] = [0.7, 0.3]
    prop['t'] = 273.15 + 35
    prop['Q'] = 1000000000##################################################################### tp be solved #########

    #single processing
    SetMultiprocessing.off()
    print(_test(prop))

    #multiprocessing
    SetMultiprocessing.on()
    print(_test(prop))

    print ('test')
    print (SetMultiprocessing())

