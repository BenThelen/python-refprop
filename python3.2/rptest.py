#-------------------------------------------------------------------------------
#Name:              rptest
#Purpose:           test module for refprop and multiRP
#
#Author:            Thelen, B.J.
#                   thelen_ben@yahoo.com
#-------------------------------------------------------------------------------
'''Allow refprop and multiRP module functional test of all functions'''

from decimal import Decimal
import refprop, platform

def settest(test):
    '''set test module
    'refprop' or 'multiRP'
    and execute test run'''
    if test == 'refprop':
        import refprop as rp
        _maintest(rp)
    elif test == 'multiRP':
        import multiRP as rp
        _maintest(rp)

#main test def. for usage at refprop and multiRP
def _maintest(rp):
    #examples and test setup
    rp.SetErrorDebug.off() #turn on =>> for testing purpose

    if rp.test(): #if True; rptest =>>for testing purpose
        print('refprop installed correctely')

        print('test results')
        print(rp.testresult)

        print('fluidlib')
        rp.fluidlib()
        print('\n')

        prop = rp.setup('def', 'air',)
        print('setup air')
        print(prop, '\n')

        x = prop['x']

        print('critp(x)')
        print(rp.critp(x), '\n')

        print('setup water ammonia')
        print(rp.setup('def', 'water', 'ammonia',), '\n')

        #alternative setup input
        rp.setup('def', ['water', 'ammonia'],)

        x = [0.5, 0.3]
        rp.normalize(x)

        prop = rp.critp(x)
        prop = rp.therm(prop['tcrit'], prop['Dcrit'], x)
        print('therm')
        print(prop, '\n')

        p = prop['p']

        print('therm2')
        print(rp.therm2(prop['t'], prop['D'], x), '\n')

        print('therm0')
        print(rp.therm0(prop['t'], prop['D'], x), '\n')

        print('residual')
        print(rp.residual(prop['t'], prop['D'], x), '\n')

        print('entro')
        print(rp.entro(prop['t'], prop['D'], x), '\n')

        print('enthal')
        print(rp.enthal(prop['t'], prop['D'], x), '\n')

        print('ag')
        print(rp.ag(prop['t'], prop['D'], x), '\n')

        print('cvcp')
        print(rp.cvcp(prop['t'], prop['D'], x), '\n')

        print('dddp')
        print(rp.dddp(prop['t'], prop['D'], x), '\n')

        print('dddt')
        print(rp.dddt(prop['t'], prop['D'], x), '\n')

        print('dhd1')
        print(rp.dhd1(prop['t'], prop['D'], x), '\n')

        print('dpdd')
        print(rp.dpdd(prop['t'], prop['D'], x), '\n')

        print('dpdd2')
        print(rp.dpdd2(prop['t'], prop['D'], x), '\n')

        print('dpdt')
        print(rp.dpdt(prop['t'], prop['D'], x), '\n')

        D = prop['D']

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('dcdt')
            print(rp.dcdt(prop['t'], x), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('dcdt2')
            print(rp.dcdt2(prop['t'], x), '\n')

        print('fgcty')
        print(rp.fgcty(prop['t'], D, x), '\n')

        print('gibbs')
        print(rp.gibbs(prop['t'], prop['D'], x), '\n')

        #~ print('fgcty2')
        #~ print(rp.fgcty2(prop['t'], prop['D'], x), '\n')

        prop = rp.therm3(prop['t'], prop['D'], x)
        print('therm3')
        print(prop, '\n')

        D = prop['D']

        print('virb')
        print(rp.virb(prop['t'], x), '\n')

        print('virc')
        print(rp.virc(prop['t'], x), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('vird')
            print(rp.vird(prop['t'], x), '\n')

        print('virba')
        print(rp.virba(prop['t'], x), '\n')

        print('virca')
        print(rp.virca(prop['t'], x), '\n')

        print('cvcpk')
        print(rp.cvcpk(1, prop['t'], D), '\n')

        print('dbdt')
        print(rp.dbdt(prop['t'], x), '\n')

        print('dpddk')
        print(rp.dpddk(1, prop['t'], D), '\n')

        print('dpdtk')
        print(rp.dpdtk(2, prop['t'], D), '\n')

        D = 55
        t = 373

        prop = rp.press(t, D, x)
        print('press')
        print(prop, '\n')

        #p = prop['p'] disabled due to error on ver 9.105
        p = 739450.6482243013 #from previous run

        print('purefld(1)')
        prop = rp.purefld(1)
        print(prop, '\n')

        x = [1]

        resetup_test_prop_d = prop

        print('satt')
        prop = rp.satt(t, x)
        print(prop, '\n')

        print('satp')
        prop = rp.satp(prop['p'], x)
        print(prop, '\n')

        print('satd')
        print(rp.satd(prop['Dliq'], x), '\n')

        print('sath')
        print(rp.sath(47000, x, 0), '\n')

        print('sate')
        print(rp.sate(45000, x), '\n')

        print('sats')
        print(rp.sats(50, x, 0), '\n')

        print('purefld(0)')
        print(rp.purefld(0), '\n')

        x = [0.5, 0.3]
        rp.normalize(x)

        print('csatk')
        print(rp.csatk(1, t), '\n')

        print('dptsatk')
        print(rp.dptsatk(1, t), '\n')

        print('cv2pk')
        print(rp.cv2pk(2, t, D), '\n')

        print('tprho')
        print(rp.tprho(t, p, x, 2, 1, 58), '\n')

        print('flsh, tp')
        prop = rp.flsh('tp', t, p, x)
        print(prop, '\n')

        print('flsh, th')
        print(rp.flsh('tH', t, prop['h'], x, 1), '\n')

        print('flsh, tD')
        print(rp.flsh('tD', t, 30, x), '\n')

        print('info()')
        print(rp.info(), '\n')

        print('info(2)')
        print(rp.info(2), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('rmix2')
            print(rp.rmix2(x), '\n')

        print('xmass')
        prop = rp.xmass(x)
        print(prop, '\n')

        print('xmole')
        print(rp.xmole(prop['xkg']), '\n')

        print('limitx')
        print(rp.limitx(x, 'eos', t, D, p), '\n')

        print('limitk')
        print(rp.limitk('eos', 1, t, D, p), '\n')

        print('limits')
        print(rp.limits(x), '\n')

        print('flsh, ts')
        prop = rp.flsh('ts', t, 40, x)
        print(prop, '\n')

        print('flsh, te')
        print(rp.flsh('te', t, prop['e'], x), '\n')

        print('flsh, pD')
        prop = rp.flsh('Pd', p, D, x)
        print(prop, '\n')

        print('flsh, ph')
        prop = rp.flsh('ph', p, prop['h'], x)
        print(prop, '\n')

        print('flsh, ps')
        prop = rp.flsh('ps', p, prop['s'], x)
        print(prop, '\n')

        print('flsh, pe')
        prop = rp.flsh('pE', p, prop['e'], x)
        print(prop, '\n')

        print('flsh, hs')
        prop = rp.flsh('hs', prop['h'], 45, x)
        print(prop, '\n')

        print('flsh, es')
        prop = rp.flsh('es', prop['e'], prop['s'], x)
        print(prop, '\n')

        print('flsh, hs')
        prop = rp.flsh('hs', prop['h'], 45, x)
        print(prop, '\n')

        print('flsh, es')
        print(rp.flsh('es', prop['e'], prop['s'], x), '\n')

        print('flsh, Dh')
        print(rp.flsh('DH', 20, 20000, x), '\n')

        print('flsh, Ds')
        prop = rp.flsh('Ds', 20, 50, x)
        print(prop, '\n')

        print('flsh, De')
        prop = rp.flsh('DE', 20, prop['e'], x)
        print(prop, '\n')

        print('flsh, tq')
        prop = rp.flsh('tq', t, prop['q'], x)
        print(prop, '\n')

        print('flsh, pq')
        print(rp.flsh('pq', 1200, prop['q'], x), '\n')

        prop = rp.flsh('tp', 350, 1200, x)
        print('flsh, tp')
        print(prop, '\n')
        s = prop['s']
        e = prop['e']
        h = prop['h']
        D = prop['D']
        t = prop['t']
        p = prop['p']
        Dmin = 40
        Dmax = 55

        print('flsh1, liq, ph')
        print(rp.flsh1('Ph', p, h, x, 1), '\n')

        print('getphase')
        print(rp.getphase(prop), '\n')

        print('flsh1, liq, pD')
        print(rp.flsh1('PD', p, D, x), '\n')

        print('flsh1, liq, ps')
        print(rp.flsh1('Ps', p, s, x), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, liq, th')
            print(rp.flsh1('th', t, h, x, Dmin=Dmin, Dmax=Dmax), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, liq, ts')
            print(rp.flsh1('ts', t, s, x, Dmin=Dmin, Dmax=Dmax), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, liq, te')
            print(rp.flsh1('te', t, e, x, Dmin=Dmin, Dmax=Dmax), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, liq, pe')
            print(rp.flsh1('Pe', p, e, x), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, liq, hs')
            print(rp.flsh1('hs', h, s, x, Dmin=Dmin, Dmax=Dmax), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, liq, Dh')
            print(rp.flsh1('Dh', D, h, x), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, liq, Ds')
            print(rp.flsh1('Ds', D, s, x), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, liq, De')
            print(rp.flsh1('De', D, e, x), '\n')

        prop = rp.flsh('tp', 400, 100, x)
        s = prop['s']
        e = prop['e']
        h = prop['h']
        D = prop['D']
        Dmin = 0.01
        Dmax = 0.05
        t = prop['t']
        p = prop['p']

        print('flsh1, vap, ph')
        print(rp.flsh1('Ph', p, h, x, 2), '\n')

        print('getphase')
        print(rp.getphase(prop), '\n')

        print('flsh1, vap, pD')
        print(rp.flsh1('PD', p, D, x, 2), '\n')

        print('flsh1, vap, ps')
        print(rp.flsh1('Ps', p, s, x, 2), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, vap, th')
            print(rp.flsh1('th', t, h, x, Dmin=Dmin, Dmax=Dmax), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, vap, ts')
            print(rp.flsh1('ts', t, s, x, Dmin=Dmin, Dmax=Dmax), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, vap, te')
            print(rp.flsh1('te', t, e, x, Dmin=Dmin, Dmax=Dmax), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, vap, pe')
            print(rp.flsh1('Pe', p, e, x, 2), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, vap, hs')
            print(rp.flsh1('hs', h, s, x, Dmin=Dmin, Dmax=Dmax), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, vap, Dh')
            print(rp.flsh1('Dh', D, h, x), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, vap, Ds')
            print(rp.flsh1('Ds', D, s, x), '\n')

        #unsupported in Windows
        if platform.system() == 'Linux':
            print('flsh1, vap, De')
            print(rp.flsh1('De', D, e, x), '\n')

        print('cstar')
        print(rp.cstar(t, p, 8, x), '\n')

        print('fpv')
        print(rp.fpv(t, D, p, x), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('excess')
            print(rp.excess(t, p, x, kph=2), '\n')

        prop = rp.flsh('pq', 1200, 0.65, x)
        D = prop['D']
        Dliq = prop['Dliq']
        Dvap = prop['Dvap']
        xliq = prop['xliq']
        xvap = prop['xvap']
        e = prop['e']
        h = prop['h']
        s = prop['s']
        q = prop['q']
        p = prop['p']
        t = prop['t']

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('tpfl2')
            print(rp.flsh2('tp', t, p, x), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('Dhfl2')
            print(rp.flsh2('Dh', D, h, x), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('Dsfl2')
            print(rp.flsh2('Ds', D, s, x), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('Defl2')
            print(rp.flsh2('De', D, e, x), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('thfl2')
            print(rp.flsh2('th', t, h, x, ksat=0), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('tsfl2')
            print(rp.flsh2('ts', t, s, x, ksat=0), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('tefl2')
            print(rp.flsh2('te', t, e, x, ksat=0), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('tDfl2')
            print(rp.flsh2('tD', t, D, x, ksat=0), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('pDfl2')
            print(rp.flsh2('pD', p, D, x, ksat=0), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('phfl2')
            print(rp.flsh2('ph', p, h, x, ksat=0), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('psfl2')
            print(rp.flsh2('ps', p, s, x, ksat=0), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('pefl2')
            print(rp.flsh2('pe', p, e, x, ksat=0), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('tqfl2')
            print(rp.flsh2('tq', t, q, x, ksat=0), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('pqfl2')
            print(rp.flsh2('pq', p, q, x, ksat=0), '\n')

        #function not supported in Windows
        #~ if platform.system() == 'Linux':
            #~ print('Dqfl2')
            #~ print(rp.flsh2('Dq', D, q, x), '\n')

        prop = rp.flsh('tp', 340, 100, x)
        t = prop['t']
        Dliq = prop['Dliq']
        Dvap = prop['Dvap']
        xliq = prop['xliq']
        xvap = prop['xvap']

        print('qmass')
        prop = rp.qmass(prop['q'], xliq, xvap)
        print(prop, '\n')

        print('qmole')
        print(rp.qmole(prop['qkg'], prop['xlkg'], prop['xvkg']), '\n')

        print('wmol')
        print(rp.wmol(x), '\n')

        prop = rp.flsh('tp', 340, 100, x)

        print('dielec')
        print(rp.dielec(prop['t'], prop['D'], x), '\n')

        print('surten')
        print(rp.surten (t, Dliq, Dvap, xliq, xvap), '\n')

        print('surft')
        print(rp.surft(240, x), '\n')

        rp.setup('def', 'water')

        print('meltt')
        print(rp.meltt(273.15, [1]), '\n')

        print('meltp')
        print(rp.meltp(100, [1]), '\n')

        print('sublt')
        print(rp.sublt(273.15, [1]), '\n')

        print('sublp')
        print(rp.sublp(0.1, [1]), '\n')

        rp.setup('def', 'butane', 'ethane', 'propane', 'methane',)
        x = [0.5, 0.15, 0.3, 0.05]
        rp.setref('nbp')
        prop = rp.flsh('tp', 260, 200, x)
        D = prop['D']
        print('trnprp, setref NBP')
        print(rp.trnprp(260, D, x), '\n')

        print('B12')
        print(rp.b12(260, x), '\n')

        print('chempot')
        print(rp.chempot(260, D, x), '\n')

        print('fugcof')
        print(rp.fugcof(260, D, x), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('phiderv')
            print(rp.phiderv(2, 1, 260, D, x), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('getmod')
            print(rp.getmod(1, 'EOS'), '\n')

        rp.setmod('tcx', 'ecs', ['tc2', 'tc1', 'tc2', 'tc2'])
        rp.setup('def', 'butane', 'ethane', 'propane', 'methane',)
        x = [0.5, 0.15, 0.3, 0.05]
        rp.setref('nbp')
        prop = rp.flsh('tp', 260, 200, x)
        print('trnprp, setref NBP, setmod [tcx, ecs, tc2, tc1, tc2, tc2]')
        print(rp.trnprp(260, prop['D'], x), '\n')

        #function not supported in Windows
        if platform.system() == 'Linux':
            print('getmod')
            print(rp.getmod(3, 'tcx'), '\n')

        rp.setref('oth', 1, [1], 0, 0, 273, 100)
        print('setref = OTH')
        prop = rp.flsh('tp', 260, 200, x)
        print(prop, '\n')

        resetup_test_prop_a = prop

        rp.setref('???', 1, [1], 0, 0, 373, 100)
        print('setref = ???')
        prop = rp.flsh('tp', 260, 200, x)
        print(prop, '\n')

        resetup_test_prop_b = prop

        print('name')
        print(rp.name(1), '\n')

        rp.setup('def', 'butane', 'ethane', 'propane', 'methane',)
        x = [0.5, 0.15, 0.3, 0.05]
        print('getktv')
        prop = rp.getktv(1, 3)
        print(prop, '\n')

        print('setktv')
        prop = rp.setktv(1, 3, 'lin', prop['fij'], prop['hfmix'],)
        print(prop, '\n')

        resetup_test_prop_c = prop

        print('reset setktv')
        print(rp.setktv(1, 2, 'rst'), '\n')

        print('getfij')
        print(rp.getfij('LIN'), '\n')

        print('resetup_test_prop, setref, setmod')
        print(resetup_test_prop_a, '\n')

        print('resetup')
        print(rp.resetup(resetup_test_prop_a), '\n')

        print('resetup_test_prop, setref(???), setmod')
        print(resetup_test_prop_b, '\n')

        print('resetup')
        print(rp.resetup(resetup_test_prop_b), '\n')

        print('resetup_test_prop, setktv')
        print(resetup_test_prop_c, '\n')

        print('resetup')
        print(rp.resetup(resetup_test_prop_c), '\n')

        print('resetup_test_prop, purefld')
        print(resetup_test_prop_d, '\n')

        print('resetup')
        print(rp.resetup(resetup_test_prop_d), '\n')

        #normalize([0.2, 0.2, 0.1, 0.1])
        print('normalize')
        print(rp.normalize([0.2, 0.2, 0.1, 0.1]), '\n')

        #setup_details
        print('setup_details')
        print(rp.setup_details({'hfld': ['BUTANE', 'ETHANE', 'PROPANE', 'METHANE'],
                                    'D': 0.21683907260570098,
                                    'Dvap': 0.09664613429889905, 'hfmix': 'HMX.BNC',
                                    'setmod': {'hcomp': ['TC2', 'TC1', 'TC2', 'TC2'],
                                                  'htype': 'TCX', 'hmix': 'ECS'},
                                    'cp': -9999980.0,
                                    'xliq': [Decimal('0.7125650648765283717349528049'),
                                                Decimal('0.04065955068790887177072495080'),
                                                Decimal('0.2449672538076863186375885862'),
                                                Decimal('0.001808130627876437856733658079')],
                                    'xvap': [Decimal('0.2304027911956556081031262882'),
                                                Decimal('0.2886769748808782463382744488'),
                                                Decimal('0.3697982730402927396744896960'),
                                                Decimal('0.1111219608831734058841095670')],
                                    'x': [0.5, 0.15, 0.3, 0.05], 'e': -13828.39837781548,
                                    'h': -12906.055381248256, 'nc': 4,
                                    'Dliq': 11.150114864150222, 'cv': -9999980.0,
                                    'q': 0.4408579356823604, 'p': 200.0,
                                    's': -44.047682476988044, 't': 260.0, 'w': -9999980.0,
                                    'kph': 1, 'setref': {'p0': 100, 'ixflag': 1, 'h0': 0,
                                                                's0': 0, 't0': 273,
                                                                'hrf': ['OTH', '???']},
                                    'hrf': 'DEF'}), '\n')

        #gerg04
        print('gerg04 = 1')
        rp.gerg04(1)
        print(rp.setup('def', 'butane', 'ethane', 'propane'), '\n')

        #reset gerg04
        print('gerg04 = 0')
        rp.gerg04(0)
        print(rp.setup('def', 'butane', 'ethane', 'propane'), '\n')

        #preos
        print('preos = 2')
        print(rp.preos(2), '\n')

        print('preos = -1')
        print(rp.preos(-1), '\n')

        print('preos = 0')
        print(rp.preos(0), '\n')

        print('preos = -1')
        print(rp.preos(-1), '\n')

        #setaga
        print('setaga')
        print(rp.setaga(), '\n')

        #unsetaga
        print('unsetaga')
        print(rp.unsetaga(), '\n')

        #setup_settings
        print('setup_setting')
        print(rp.setup_setting(), '\n')
