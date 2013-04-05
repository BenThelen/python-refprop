#-------------------------------------------------------------------------------
#Name:              rptest
#Purpose:           test module for refprop and multiRP
#
#Author:            Thelen, B.J.
#                   thelen_ben@yahoo.com
#-------------------------------------------------------------------------------
u'''Allow refprop and multiRP module functional test of all functions'''

from decimal import Decimal
import refprop, platform

def settest(test):
    u'''set test module
    'refprop' or 'multiRP'
    and execute test run'''
    if test == u'refprop':
        import refprop as rp
        _maintest(rp)
    elif test == u'multiRP':
        import multiRP as rp
        _maintest(rp)

#main test def. for usage at refprop and multiRP
def _maintest(rp):
    #examples and test setup
    rp.SetErrorDebug.off() #turn on =>> for testing purpose

    if rp.test(): #if True; rptest =>>for testing purpose
        print u'refprop installed correctely'

        print u'test results'
        print rp.testresult

        print u'fluidlib'
        rp.fluidlib()
        print u'\n'

        prop = rp.setup(u'def', u'air',)
        print u'setup air'
        print prop, u'\n'

        x = prop[u'x']

        print u'critp(x)'
        print rp.critp(x), u'\n'

        print u'setup water ammonia'
        print rp.setup(u'def', u'water', u'ammonia',), u'\n'

        #alternative setup input
        rp.setup(u'def', [u'water', u'ammonia'],)

        x = [0.5, 0.3]
        rp.normalize(x)

        prop = rp.critp(x)
        prop = rp.therm(prop[u'tcrit'], prop[u'Dcrit'], x)
        print u'therm'
        print prop, u'\n'

        p = prop[u'p']

        print u'therm2'
        print rp.therm2(prop[u't'], prop[u'D'], x), u'\n'

        print u'therm0'
        print rp.therm0(prop[u't'], prop[u'D'], x), u'\n'

        print u'residual'
        print rp.residual(prop[u't'], prop[u'D'], x), u'\n'

        print u'entro'
        print rp.entro(prop[u't'], prop[u'D'], x), u'\n'

        print u'enthal'
        print rp.enthal(prop[u't'], prop[u'D'], x), u'\n'

        print u'ag'
        print rp.ag(prop[u't'], prop[u'D'], x), u'\n'

        print u'cvcp'
        print rp.cvcp(prop[u't'], prop[u'D'], x), u'\n'

        print u'dddp'
        print rp.dddp(prop[u't'], prop[u'D'], x), u'\n'

        print u'dddt'
        print rp.dddt(prop[u't'], prop[u'D'], x), u'\n'

        print u'dhd1'
        print rp.dhd1(prop[u't'], prop[u'D'], x), u'\n'

        print u'dpdd'
        print rp.dpdd(prop[u't'], prop[u'D'], x), u'\n'

        print u'dpdd2'
        print rp.dpdd2(prop[u't'], prop[u'D'], x), u'\n'

        print u'dpdt'
        print rp.dpdt(prop[u't'], prop[u'D'], x), u'\n'

        D = prop[u'D']

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'dcdt'
            print rp.dcdt(prop[u't'], x), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'dcdt2'
            print rp.dcdt2(prop[u't'], x), u'\n'

        print u'fgcty'
        print rp.fgcty(prop[u't'], D, x), u'\n'

        print u'gibbs'
        print rp.gibbs(prop[u't'], prop[u'D'], x), u'\n'

        #~ print('fgcty2')
        #~ print(rp.fgcty2(prop['t'], prop['D'], x), '\n')

        prop = rp.therm3(prop[u't'], prop[u'D'], x)
        print u'therm3'
        print prop, u'\n'

        D = prop[u'D']

        print u'virb'
        print rp.virb(prop[u't'], x), u'\n'

        print u'virc'
        print rp.virc(prop[u't'], x), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'vird'
            print rp.vird(prop[u't'], x), u'\n'

        print u'virba'
        print rp.virba(prop[u't'], x), u'\n'

        print u'virca'
        print rp.virca(prop[u't'], x), u'\n'

        print u'cvcpk'
        print rp.cvcpk(1, prop[u't'], D), u'\n'

        print u'dbdt'
        print rp.dbdt(prop[u't'], x), u'\n'

        print u'dpddk'
        print rp.dpddk(1, prop[u't'], D), u'\n'

        print u'dpdtk'
        print rp.dpdtk(2, prop[u't'], D), u'\n'

        D = 55
        t = 373

        prop = rp.press(t, D, x)
        print u'press'
        print prop, u'\n'

        #p = prop['p'] disabled due to error on ver 9.105
        p = 739450.6482243013 #from previous run

        print u'purefld(1)'
        prop = rp.purefld(1)
        print prop, u'\n'

        x = [1]

        resetup_test_prop_d = prop

        print u'satt'
        prop = rp.satt(t, x)
        print prop, u'\n'

        print u'satp'
        prop = rp.satp(prop[u'p'], x)
        print prop, u'\n'

        print u'satd'
        print rp.satd(prop[u'Dliq'], x), u'\n'

        print u'sath'
        print rp.sath(47000, x, 0), u'\n'

        print u'sate'
        print rp.sate(45000, x), u'\n'

        print u'sats'
        print rp.sats(50, x, 0), u'\n'

        print u'purefld(0)'
        print rp.purefld(0), u'\n'

        x = [0.5, 0.3]
        rp.normalize(x)

        print u'csatk'
        print rp.csatk(1, t), u'\n'

        print u'dptsatk'
        print rp.dptsatk(1, t), u'\n'

        print u'cv2pk'
        print rp.cv2pk(2, t, D), u'\n'

        print u'tprho'
        print rp.tprho(t, p, x, 2, 1, 58), u'\n'

        print u'flsh, tp'
        prop = rp.flsh(u'tp', t, p, x)
        print prop, u'\n'

        print u'flsh, th'
        print rp.flsh(u'tH', t, prop[u'h'], x, 1), u'\n'

        print u'flsh, tD'
        print rp.flsh(u'tD', t, 30, x), u'\n'

        print u'info()'
        print rp.info(), u'\n'

        print u'info(2)'
        print rp.info(2), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'rmix2'
            print rp.rmix2(x), u'\n'

        print u'xmass'
        prop = rp.xmass(x)
        print prop, u'\n'

        print u'xmole'
        print rp.xmole(prop[u'xkg']), u'\n'

        print u'limitx'
        print rp.limitx(x, u'eos', t, D, p), u'\n'

        print u'limitk'
        print rp.limitk(u'eos', 1, t, D, p), u'\n'

        print u'limits'
        print rp.limits(x), u'\n'

        print u'flsh, ts'
        prop = rp.flsh(u'ts', t, 40, x)
        print prop, u'\n'

        print u'flsh, te'
        print rp.flsh(u'te', t, prop[u'e'], x), u'\n'

        print u'flsh, pD'
        prop = rp.flsh(u'Pd', p, D, x)
        print prop, u'\n'

        print u'flsh, ph'
        prop = rp.flsh(u'ph', p, prop[u'h'], x)
        print prop, u'\n'

        print u'flsh, ps'
        prop = rp.flsh(u'ps', p, prop[u's'], x)
        print prop, u'\n'

        print u'flsh, pe'
        prop = rp.flsh(u'pE', p, prop[u'e'], x)
        print prop, u'\n'

        print u'flsh, hs'
        prop = rp.flsh(u'hs', prop[u'h'], 45, x)
        print prop, u'\n'

        print u'flsh, es'
        prop = rp.flsh(u'es', prop[u'e'], prop[u's'], x)
        print prop, u'\n'

        print u'flsh, hs'
        prop = rp.flsh(u'hs', prop[u'h'], 45, x)
        print prop, u'\n'

        print u'flsh, es'
        print rp.flsh(u'es', prop[u'e'], prop[u's'], x), u'\n'

        print u'flsh, Dh'
        print rp.flsh(u'DH', 20, 20000, x), u'\n'

        print u'flsh, Ds'
        prop = rp.flsh(u'Ds', 20, 50, x)
        print prop, u'\n'

        print u'flsh, De'
        prop = rp.flsh(u'DE', 20, prop[u'e'], x)
        print prop, u'\n'

        print u'flsh, tq'
        prop = rp.flsh(u'tq', t, prop[u'q'], x)
        print prop, u'\n'

        print u'flsh, pq'
        print rp.flsh(u'pq', 1200, prop[u'q'], x), u'\n'

        prop = rp.flsh(u'tp', 350, 1200, x)
        print u'flsh, tp'
        print prop, u'\n'
        s = prop[u's']
        e = prop[u'e']
        h = prop[u'h']
        D = prop[u'D']
        t = prop[u't']
        p = prop[u'p']
        Dmin = 40
        Dmax = 55

        print u'flsh1, liq, ph'
        print rp.flsh1(u'Ph', p, h, x, 1), u'\n'

        print u'getphase'
        print rp.getphase(prop), u'\n'

        print u'flsh1, liq, pD'
        print rp.flsh1(u'PD', p, D, x), u'\n'

        print u'flsh1, liq, ps'
        print rp.flsh1(u'Ps', p, s, x), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, liq, th'
            print rp.flsh1(u'th', t, h, x, Dmin=Dmin, Dmax=Dmax), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, liq, ts'
            print rp.flsh1(u'ts', t, s, x, Dmin=Dmin, Dmax=Dmax), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, liq, te'
            print rp.flsh1(u'te', t, e, x, Dmin=Dmin, Dmax=Dmax), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, liq, pe'
            print rp.flsh1(u'Pe', p, e, x), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, liq, hs'
            print rp.flsh1(u'hs', h, s, x, Dmin=Dmin, Dmax=Dmax), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, liq, Dh'
            print rp.flsh1(u'Dh', D, h, x), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, liq, Ds'
            print rp.flsh1(u'Ds', D, s, x), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, liq, De'
            print rp.flsh1(u'De', D, e, x), u'\n'

        prop = rp.flsh(u'tp', 400, 100, x)
        s = prop[u's']
        e = prop[u'e']
        h = prop[u'h']
        D = prop[u'D']
        Dmin = 0.01
        Dmax = 0.05
        t = prop[u't']
        p = prop[u'p']

        print u'flsh1, vap, ph'
        print rp.flsh1(u'Ph', p, h, x, 2), u'\n'

        print u'getphase'
        print rp.getphase(prop), u'\n'

        print u'flsh1, vap, pD'
        print rp.flsh1(u'PD', p, D, x, 2), u'\n'

        print u'flsh1, vap, ps'
        print rp.flsh1(u'Ps', p, s, x, 2), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, vap, th'
            print rp.flsh1(u'th', t, h, x, Dmin=Dmin, Dmax=Dmax), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, vap, ts'
            print rp.flsh1(u'ts', t, s, x, Dmin=Dmin, Dmax=Dmax), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, vap, te'
            print rp.flsh1(u'te', t, e, x, Dmin=Dmin, Dmax=Dmax), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, vap, pe'
            print rp.flsh1(u'Pe', p, e, x, 2), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, vap, hs'
            print rp.flsh1(u'hs', h, s, x, Dmin=Dmin, Dmax=Dmax), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, vap, Dh'
            print rp.flsh1(u'Dh', D, h, x), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, vap, Ds'
            print rp.flsh1(u'Ds', D, s, x), u'\n'

        #unsupported in Windows
        if platform.system() == u'Linux':
            print u'flsh1, vap, De'
            print rp.flsh1(u'De', D, e, x), u'\n'

        print u'cstar'
        print rp.cstar(t, p, 8, x), u'\n'

        print u'fpv'
        print rp.fpv(t, D, p, x), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'excess'
            print rp.excess(t, p, x, kph=2), u'\n'

        prop = rp.flsh(u'pq', 1200, 0.65, x)
        D = prop[u'D']
        Dliq = prop[u'Dliq']
        Dvap = prop[u'Dvap']
        xliq = prop[u'xliq']
        xvap = prop[u'xvap']
        e = prop[u'e']
        h = prop[u'h']
        s = prop[u's']
        q = prop[u'q']
        p = prop[u'p']
        t = prop[u't']

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'tpfl2'
            print rp.flsh2(u'tp', t, p, x), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'Dhfl2'
            print rp.flsh2(u'Dh', D, h, x), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'Dsfl2'
            print rp.flsh2(u'Ds', D, s, x), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'Defl2'
            print rp.flsh2(u'De', D, e, x), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'thfl2'
            print rp.flsh2(u'th', t, h, x, ksat=0), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'tsfl2'
            print rp.flsh2(u'ts', t, s, x, ksat=0), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'tefl2'
            print rp.flsh2(u'te', t, e, x, ksat=0), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'tDfl2'
            print rp.flsh2(u'tD', t, D, x, ksat=0), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'pDfl2'
            print rp.flsh2(u'pD', p, D, x, ksat=0), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'phfl2'
            print rp.flsh2(u'ph', p, h, x, ksat=0), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'psfl2'
            print rp.flsh2(u'ps', p, s, x, ksat=0), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'pefl2'
            print rp.flsh2(u'pe', p, e, x, ksat=0), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'tqfl2'
            print rp.flsh2(u'tq', t, q, x, ksat=0), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'pqfl2'
            print rp.flsh2(u'pq', p, q, x, ksat=0), u'\n'

        #function not supported in Windows
        #~ if platform.system() == 'Linux':
            #~ print('Dqfl2')
            #~ print(rp.flsh2('Dq', D, q, x), '\n')

        prop = rp.flsh(u'tp', 340, 100, x)
        t = prop[u't']
        Dliq = prop[u'Dliq']
        Dvap = prop[u'Dvap']
        xliq = prop[u'xliq']
        xvap = prop[u'xvap']

        print u'qmass'
        prop = rp.qmass(prop[u'q'], xliq, xvap)
        print prop, u'\n'

        print u'qmole'
        print rp.qmole(prop[u'qkg'], prop[u'xlkg'], prop[u'xvkg']), u'\n'

        print u'wmol'
        print rp.wmol(x), u'\n'

        prop = rp.flsh(u'tp', 340, 100, x)

        print u'dielec'
        print rp.dielec(prop[u't'], prop[u'D'], x), u'\n'

        print u'surten'
        print rp.surten (t, Dliq, Dvap, xliq, xvap), u'\n'

        print u'surft'
        print rp.surft(240, x), u'\n'

        rp.setup(u'def', u'water')

        print u'meltt'
        print rp.meltt(273.15, [1]), u'\n'

        print u'meltp'
        print rp.meltp(100, [1]), u'\n'

        print u'sublt'
        print rp.sublt(273.15, [1]), u'\n'

        print u'sublp'
        print rp.sublp(0.1, [1]), u'\n'

        rp.setup(u'def', u'butane', u'ethane', u'propane', u'methane',)
        x = [0.5, 0.15, 0.3, 0.05]
        rp.setref(u'nbp')
        prop = rp.flsh(u'tp', 260, 200, x)
        D = prop[u'D']
        print u'trnprp, setref NBP'
        print rp.trnprp(260, D, x), u'\n'

        print u'B12'
        print rp.b12(260, x), u'\n'

        print u'chempot'
        print rp.chempot(260, D, x), u'\n'

        print u'fugcof'
        print rp.fugcof(260, D, x), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'phiderv'
            print rp.phiderv(2, 1, 260, D, x), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'getmod'
            print rp.getmod(1, u'EOS'), u'\n'

        rp.setmod(u'tcx', u'ecs', [u'tc2', u'tc1', u'tc2', u'tc2'])
        rp.setup(u'def', u'butane', u'ethane', u'propane', u'methane',)
        x = [0.5, 0.15, 0.3, 0.05]
        rp.setref(u'nbp')
        prop = rp.flsh(u'tp', 260, 200, x)
        print u'trnprp, setref NBP, setmod [tcx, ecs, tc2, tc1, tc2, tc2]'
        print rp.trnprp(260, prop[u'D'], x), u'\n'

        #function not supported in Windows
        if platform.system() == u'Linux':
            print u'getmod'
            print rp.getmod(3, u'tcx'), u'\n'

        rp.setref(u'oth', 1, [1], 0, 0, 273, 100)
        print u'setref = OTH'
        prop = rp.flsh(u'tp', 260, 200, x)
        print prop, u'\n'

        resetup_test_prop_a = prop

        rp.setref(u'???', 1, [1], 0, 0, 373, 100)
        print u'setref = ???'
        prop = rp.flsh(u'tp', 260, 200, x)
        print prop, u'\n'

        resetup_test_prop_b = prop

        print u'name'
        print rp.name(1), u'\n'

        rp.setup(u'def', u'butane', u'ethane', u'propane', u'methane',)
        x = [0.5, 0.15, 0.3, 0.05]
        print u'getktv'
        prop = rp.getktv(1, 3)
        print prop, u'\n'

        print u'setktv'
        prop = rp.setktv(1, 3, u'lin', prop[u'fij'], prop[u'hfmix'],)
        print prop, u'\n'

        resetup_test_prop_c = prop

        print u'reset setktv'
        print rp.setktv(1, 2, u'rst'), u'\n'

        print u'getfij'
        print rp.getfij(u'LIN'), u'\n'

        print u'resetup_test_prop, setref, setmod'
        print resetup_test_prop_a, u'\n'

        print u'resetup'
        print rp.resetup(resetup_test_prop_a), u'\n'

        print u'resetup_test_prop, setref(???), setmod'
        print resetup_test_prop_b, u'\n'

        print u'resetup'
        print rp.resetup(resetup_test_prop_b), u'\n'

        print u'resetup_test_prop, setktv'
        print resetup_test_prop_c, u'\n'

        print u'resetup'
        print rp.resetup(resetup_test_prop_c), u'\n'

        print u'resetup_test_prop, purefld'
        print resetup_test_prop_d, u'\n'

        print u'resetup'
        print rp.resetup(resetup_test_prop_d), u'\n'

        #normalize([0.2, 0.2, 0.1, 0.1])
        print u'normalize'
        print rp.normalize([0.2, 0.2, 0.1, 0.1]), u'\n'

        #setup_details
        print u'setup_details'
        print rp.setup_details({u'hfld': [u'BUTANE', u'ETHANE', u'PROPANE', u'METHANE'],
                                    u'D': 0.21683907260570098,
                                    u'Dvap': 0.09664613429889905, u'hfmix': u'HMX.BNC',
                                    u'setmod': {u'hcomp': [u'TC2', u'TC1', u'TC2', u'TC2'],
                                                  u'htype': u'TCX', u'hmix': u'ECS'},
                                    u'cp': -9999980.0,
                                    u'xliq': [Decimal(u'0.7125650648765283717349528049'),
                                                Decimal(u'0.04065955068790887177072495080'),
                                                Decimal(u'0.2449672538076863186375885862'),
                                                Decimal(u'0.001808130627876437856733658079')],
                                    u'xvap': [Decimal(u'0.2304027911956556081031262882'),
                                                Decimal(u'0.2886769748808782463382744488'),
                                                Decimal(u'0.3697982730402927396744896960'),
                                                Decimal(u'0.1111219608831734058841095670')],
                                    u'x': [0.5, 0.15, 0.3, 0.05], u'e': -13828.39837781548,
                                    u'h': -12906.055381248256, u'nc': 4,
                                    u'Dliq': 11.150114864150222, u'cv': -9999980.0,
                                    u'q': 0.4408579356823604, u'p': 200.0,
                                    u's': -44.047682476988044, u't': 260.0, u'w': -9999980.0,
                                    u'kph': 1, u'setref': {u'p0': 100, u'ixflag': 1, u'h0': 0,
                                                                u's0': 0, u't0': 273,
                                                                u'hrf': [u'OTH', u'???']},
                                    u'hrf': u'DEF'}), u'\n'

        #gerg04
        print u'gerg04 = 1'
        rp.gerg04(1)
        print rp.setup(u'def', u'butane', u'ethane', u'propane'), u'\n'

        #reset gerg04
        print u'gerg04 = 0'
        rp.gerg04(0)
        print rp.setup(u'def', u'butane', u'ethane', u'propane'), u'\n'

        #preos
        print u'preos = 2'
        print rp.preos(2), u'\n'

        print u'preos = -1'
        print rp.preos(-1), u'\n'

        print u'preos = 0'
        print rp.preos(0), u'\n'

        print u'preos = -1'
        print rp.preos(-1), u'\n'

        #setaga
        print u'setaga'
        print rp.setaga(), u'\n'

        #unsetaga
        print u'unsetaga'
        print rp.unsetaga(), u'\n'

        #setup_settings
        print u'setup_setting'
        print rp.setup_setting(), u'\n'
