#!/usr/bin/env python
from lib.nekTestCase import *
from unittest import skip

import re

####################################################################
#  annulus_2d; cyl2.rea
####################################################################

class Annulus2d(NekTestCase):
    example_subdir  = 'annulus_2d'
    case_name        = 'cyl2'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '10',
            lxd       = '15',
            lx2       = 'lx1-2',
            lelg      = '100',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
            lfdm      = '0',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek()

        gmres = self.get_value_from_log('gmres', column=-7,row=1)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=20., label='gmres')

        vmax = self.get_value_from_log('tmax', column=-3, row=-1)
        self.assertAlmostEqualDelayed(vmax, target_val=2.369E-01, delta=5E-04, label='vmax')
        
        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek()

        gmres = self.get_value_from_log('gmres', column=-6,row=1)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=20., label='gmres')

        vmax = self.get_value_from_log('tmax', column=-3, row=-1)
        self.assertAlmostEqualDelayed(vmax, target_val=2.369E-01, delta=5E-04, label='vmax')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  benard: ray_9.rea
####################################################################

class Benard_Ray9(NekTestCase):
    example_subdir = 'benard'
    case_name = 'ray_9'
 
    def setUp(self):
        self.size_params = dict (
            ldim      = '2',
            lx1       = '12',
            lxd       = '16',
            lx2       = 'lx1-2',
            lelg      = '500',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
 
        self.build_tools(['clean','genmap'])
        self.run_genmap(rea_file='ray_9')
 
    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek(usr_file='ray_9')
        self.run_nek(rea_file='ray_9')
 
        gmres = self.get_value_from_log('gmres ', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=17., label='gmres')
 
        rac = self.get_value_from_log('converged_rac', column=-2)
        self.assertAlmostEqualDelayed(rac, target_val=1707.79, delta=0.01, label='rac')
 
    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek(usr_file='ray_9')
        self.run_nek(rea_file='ray_9')
 
        gmres = self.get_value_from_log('gmres ', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=11., label='gmres')
 
        rac = self.get_value_from_log('converged_rac', column=-2)
        self.assertAlmostEqualDelayed(rac, target_val=1707.79, delta=0.01, label='rac')
 
    def tearDown(self):
        self.move_logs()

####################################################################
#  blasius: blasius.rea
####################################################################

class Blasius(NekTestCase):
    example_subdir  = 'blasius'
    case_name        = 'blasius'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '1000',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=162., label='gmres')

        delta = self.get_value_from_log('delta', column=-5, row=-1)
        self.assertAlmostEqualDelayed(delta, target_val=1.26104, delta=1e-05, label='delta')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=125., label='gmres')

        delta = self.get_value_from_log('delta', column=-5, row=-1)
        self.assertAlmostEqualDelayed(delta, target_val=1.26104, delta=1e-05, label='delta')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  CMT/inviscid_vortex: pvort.rea
####################################################################

class CmtInviscidVortex(NekTestCase):
    example_subdir = os.path.join('CMT', 'inv_vort')
    case_name = 'pvort'

    def diff_l2norms(self):
        def get_line(filename, line_num=0):
            with open(filename) as f:
                line = f.readlines()[line_num]
            return [float(x) for x in line.split()[1:]]

        cls = self.__class__
        test_vals = get_line(os.path.join(self.examples_root, cls.example_subdir, 'l2norms.dat'))
        ref_vals = get_line(os.path.join(self.examples_root, cls.example_subdir, 'l2norms.dat.ref'))
        for t, r in zip(test_vals, ref_vals):
            self.assertAlmostEqual(t, r, delta=0.1*r,
                msg='FAILURE: Last line of l2norms.dat differed from reference values by > 10%\n  test vals:{0}\n  ref vals: {1}'.format(test_vals, ref_vals))
        print('SUCCESS: Last line of l2norms.dat was within 10% of reference values\n  test vals:{0}\n  ref vals: {1}'.format(test_vals, ref_vals))

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '25',
            lxd       = '36',
            lx2       = 'lx1-0',
            lelg      = '50',
            ldimt     = '3',
            toteq     = '5',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.config_size()
        self.build_tools(['clean','genmap'])
        self.run_genmap()

        cls = self.__class__
        try:
            os.remove(os.path.join(self.examples_root, cls.example_subdir, 'l2norms.dat'))
        except OSError:
            pass

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.build_nek(opts={'PPLIST':'CMTNEK'})
        self.run_nek(step_limit=1000)
        self.diff_l2norms()

    def tearDown(self):
        self.move_logs()

# ####################################################################
# #  cone/{cone016, cone064, cone256}: cone.rea
# ####################################################################

class Cone16(NekTestCase):
    example_subdir = os.path.join('cone','cone016')
    case_name = 'cone'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '9',
            lxd       = '14',
            lx2       = 'lx1-2',
            lelg      = '100',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genbox','genmap'])
        self.run_genbox(box_file='cone016')
        self.run_genmap(rea_file='box',tol='0.01')
        self.mvn('box', 'cone')
        
    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek()

        umin = self.get_value_from_log('Tmax', column=-3)
        self.assertAlmostEqualDelayed(umin, target_val=-2.3032E-02, delta=5E-05, label='Umin')
        
        umax = self.get_value_from_log('Tmax', column=-2)
        self.assertAlmostEqualDelayed(umax, target_val=8.5065E-01, delta=5E-05, label='Umax')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek()

        umin = self.get_value_from_log('Tmax', column=-3)
        self.assertAlmostEqualDelayed(umin, target_val=-2.3032E-02, delta=5E-05, label='Umin')

        umax = self.get_value_from_log('Tmax', column=-2)
        self.assertAlmostEqualDelayed(umax, target_val=8.5065E-01, delta=5E-05, label='Umax')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()    
        
class Cone64(NekTestCase):
    example_subdir = os.path.join('cone','cone064')
    case_name = 'cone'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '5',
            lxd       = '8',
            lx2       = 'lx1-2',
            lelg      = '100',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genbox','genmap'])
        self.run_genbox(box_file='cone064')
        self.run_genmap(rea_file='box',tol='0.01')
        self.mvn('box', 'cone')
        
    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek()

        umin = self.get_value_from_log('Tmax', column=-3)
        self.assertAlmostEqualDelayed(umin, target_val=-1.2663E-01, delta=5E-04, label='Umin')

        umax = self.get_value_from_log('Tmax', column=-2)
        self.assertAlmostEqualDelayed(umax, target_val=7.9285E-01, delta=5E-04, label='Umax')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek()

        umin = self.get_value_from_log('Tmax', column=-3)
        self.assertAlmostEqualDelayed(umin, target_val=-1.2663E-01, delta=5E-04, label='Umin')

        umax = self.get_value_from_log('Tmax', column=-2)
        self.assertAlmostEqualDelayed(umax, target_val=7.9285E-01, delta=5E-04, label='Umax')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class Cone256(NekTestCase):
    example_subdir = os.path.join('cone','cone256')
    case_name = 'cone'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '3',
            lxd       = '5',
            lx2       = 'lx1-2',
            lelg      = '300',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genbox','genmap'])
        self.run_genbox(box_file='cone256')
        self.run_genmap(rea_file='box',tol='0.01')
        self.mvn('box', 'cone')
        
    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek()

        umin = self.get_value_from_log('Tmax', column=-3)
        self.assertAlmostEqualDelayed(umin, target_val=-1.6392E-01, delta=5E-04, label='Umin')
        
        umax = self.get_value_from_log('Tmax', column=-2)
        self.assertAlmostEqualDelayed(umax, target_val=7.4924E-01, delta=5E-04, label='Umax')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek()

        umin = self.get_value_from_log('Tmax', column=-3)
        self.assertAlmostEqualDelayed(umin, target_val=-1.6392E-01, delta=5E-04, label='Umin')
        
        umax = self.get_value_from_log('Tmax', column=-2)
        self.assertAlmostEqualDelayed(umax, target_val=7.4924E-01, delta=5E-04, label='Umax')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()
        
####################################################################
#  conj_ht: conj_ht.rea
####################################################################

class ConjHt(NekTestCase):
    example_subdir  = 'conj_ht'
    case_name        = 'conj_ht'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '4',
            lxd       = '7',
            lx2       = 'lx1',
            lelg      = '100',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres', column=-7,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=46., label='gmres')

        tmax = self.get_value_from_log('tmax', column=-2, row=-1)
        self.assertAlmostEqualDelayed(tmax, target_val=13.1073, delta=1E-06, label='tmax')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=26., label='gmres')

        tmax = self.get_value_from_log('tmax', column=-2, row=-1)
        self.assertAlmostEqualDelayed(tmax, target_val=1.31190E+01, delta=1E-06, label='tmax')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  cyl_restart: ca.rea, cb.rea, pa.rea, pb.rea
####################################################################

class CylRestart_C(NekTestCase):
    example_subdir  = 'cyl_restart'
    case_name        = 'ca'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '6',
            lxd       = '8',
            lx2       = 'lx1-2',
            lelg      = '200',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')
        self.run_genmap(rea_file='cb', tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        import os.path
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres_a = self.get_value_from_log('gmres', column=-7, row=-1,)
        ref_val_x = self.get_value_from_log('1dragx', column=-4, row=-1)
        ref_val_y = self.get_value_from_log('1dragy', column=-4, row=-1)

        self.build_nek(usr_file='cb')
        self.run_nek(rea_file='cb',step_limit=None)
        
        logfl = os.path.join(
                self.examples_root,
                self.example_subdir,
                '{0}.log.{1}{2}'.format('cb', self.mpi_procs, self.log_suffix)
            )

        gmres_b = self.get_value_from_log('gmres', column=-7, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(gmres_b, target_val=gmres_a, delta=1., label='gmres b')

        test_val_x = self.get_value_from_log('1dragx', column=-4, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(test_val_x, target_val=ref_val_x, delta=1E-06, label='1dragx')

        test_val_y = self.get_value_from_log('1dragy', column=-4, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(test_val_y, target_val=ref_val_y, delta=1E-06, label='1dragy')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        import os.path
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres_a = self.get_value_from_log('gmres', column=-6, row=-1,)
        ref_val_x = self.get_value_from_log('1dragx', column=-4, row=-1)
        ref_val_y = self.get_value_from_log('1dragy', column=-4, row=-1)

        self.build_nek(usr_file='cb')
        self.run_nek(rea_file='cb',step_limit=None)
        
        logfl = os.path.join(
                self.examples_root,
                self.example_subdir,
                '{0}.log.{1}{2}'.format('cb', self.mpi_procs, self.log_suffix)
            )

        gmres_b = self.get_value_from_log('gmres', column=-6, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(gmres_b, target_val=gmres_a, delta=1., label='gmres b')

        test_val_x = self.get_value_from_log('1dragx', column=-4, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(test_val_x, target_val=ref_val_x, delta=1E-06, label='1dragx')

        test_val_y = self.get_value_from_log('1dragy', column=-4, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(test_val_y, target_val=ref_val_y, delta=1E-06, label='1dragy')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class CylRestart_P(NekTestCase):
    example_subdir  = 'cyl_restart'
    case_name        = 'pa'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '6',
            lxd       = '8',
            lx2       = 'lx1-2',
            lelg      = '200',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')
        self.run_genmap(rea_file='pb', tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        import os.path
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres_a = self.get_value_from_log('gmres', column=-7, row=-1,)
        ref_val_x = self.get_value_from_log('1dragx', column=-4, row=-1)
        ref_val_y = self.get_value_from_log('1dragy', column=-4, row=-1)

        self.build_nek(usr_file='pb')
        self.run_nek(rea_file='pb',step_limit=None)
        
        logfl = os.path.join(
                self.examples_root,
                self.example_subdir,
                '{0}.log.{1}{2}'.format('pb', self.mpi_procs, self.log_suffix)
            )

        gmres_b = self.get_value_from_log('gmres', column=-7, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(gmres_b, target_val=gmres_a, delta=5., label='gmres b')

        test_val_x = self.get_value_from_log('1dragx', column=-4, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(test_val_x, target_val=ref_val_x, delta=1E-06, label='1dragx')

        test_val_y = self.get_value_from_log('1dragy', column=-4, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(test_val_y, target_val=ref_val_y, delta=1E-06, label='1dragy')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        import os.path
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres_a = self.get_value_from_log('gmres', column=-6, row=-1,)
        ref_val_x = self.get_value_from_log('1dragx', column=-4, row=-1)
        ref_val_y = self.get_value_from_log('1dragy', column=-4, row=-1)

        self.build_nek(usr_file='pb')
        self.run_nek(rea_file='pb',step_limit=None)
        
        logfl = os.path.join(
                self.examples_root,
                self.example_subdir,
                '{0}.log.{1}{2}'.format('pb', self.mpi_procs, self.log_suffix)
            )

        gmres_b = self.get_value_from_log('gmres', column=-6, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(gmres_b, target_val=gmres_a, delta=5., label='gmres b')

        test_val_x = self.get_value_from_log('1dragx', column=-4, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(test_val_x, target_val=ref_val_x, delta=1E-06, label='1dragx')

        test_val_y = self.get_value_from_log('1dragy', column=-4, row=-1, logfile=logfl)
        self.assertAlmostEqualDelayed(test_val_y, target_val=ref_val_y, delta=1E-06, label='1dragy')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

#####################################################################
#   eddy_neknek: eddy_neknek.rea
#####################################################################

class Eddy_Neknek(NekTestCase):
    example_subdir  = 'eddy_neknek'
    case_name       = 'eddy_uv'
 
    def setUp(self):
 
        self.size_params = dict(
            ldim='2',
            lx1='8',
            lxd='12',
            lx2='lx1-2',
            lelg='1000',
            lpert='1',
            nsessmax='2',
        )
 
        self.build_tools(['clean','genmap'])
 
    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        from lib.nekBinRun import run_neknek
        from re import sub
 
        cls = self.__class__
        cwd = os.path.join(self.examples_root, cls.example_subdir)
 
        # Tweak the .rea files and run genmap
        for rea_file in ('inside', 'outside'):
            rea_path = os.path.join(cwd, rea_file + '.rea')
            with open(rea_path, 'r') as f:
                lines = [sub(r'^.*DIVERGENCE$', '      1.0000000E-06     p21 DIVERGENCE', l) for l in f]
            with open(rea_path, 'w') as f:
                f.writelines(lines)
            self.run_genmap(os.path.join(cwd, rea_file),tol='0.2')
 
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek(opts={'PPLIST':'NEKNEK'})
        run_neknek(
            cwd = cwd,
            inside = 'inside',
            outside = 'outside',
            np_inside = 1,
            np_outside = 1,
            step_limit = 1000,
            log_suffix = self.log_suffix,
            verbose = self.verbose,
        )
 
        logfile  = os.path.join(cwd, '{inside}{np_in}.{outside}{np_out}.log{sfx}'.format(
            inside = 'inside',
            outside = 'outside',
            np_in = 1,
            np_out = 1,
            sfx = self.log_suffix
        ))
 
        xerr_inside = self.get_value_from_log('X err  inside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_inside, target_val=7.163001E-04, delta=1E-05, label='X err  inside')
 
        xerr_global = self.get_value_from_log('X err   global', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_global, target_val=8.580050E-04, delta=1E-05, label='X err   global')
 
        xerr_outside = self.get_value_from_log('X err  outside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_outside, target_val=8.580050E-04, delta=1E-05, label='X err  outside')
 
        yerr_inside = self.get_value_from_log('Y err  inside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_inside, target_val=9.012947E-04, delta=1E-05, label='Y err  inside')
 
        yerr_global = self.get_value_from_log('Y err   global', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_global, target_val=9.877146E-04, delta=1E-05, label='Y err   global')
 
        yerr_outside = self.get_value_from_log('Y err  outside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_outside, target_val=9.877146E-04, delta=1E-05, label='Y err  outside')
 
        self.assertDelayedFailures()
 
    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        from lib.nekBinRun import run_neknek
        from re import sub
 
        cls = self.__class__
        cwd = os.path.join(self.examples_root, cls.example_subdir)
 
        # Tweak the .rea files and run genmap
        for rea_file in ('inside', 'outside'):
            rea_path = os.path.join(cwd, rea_file + '.rea')
            with open(rea_path, 'r') as f:
                lines = [sub(r'^.*DIVERGENCE$', '      1.0000000E-11     p21 DIVERGENCE', l) for l in f]
            with open(rea_path, 'w') as f:
                f.writelines(lines)
            self.run_genmap(os.path.join(cwd, rea_file),tol='0.2')
 
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek(opts={'PPLIST':'NEKNEK'})
        run_neknek(
            cwd = cwd,
            inside = 'inside',
            outside = 'outside',
            np_inside = 1,
            np_outside = 1,
            step_limit = 1000,
            log_suffix = self.log_suffix,
            verbose = self.verbose,
        )
 
        logfile  = os.path.join(cwd, '{inside}{np_in}.{outside}{np_out}.log{sfx}'.format(
            inside = 'inside',
            outside = 'outside',
            np_in = 1,
            np_out = 1,
            sfx = self.log_suffix
        ))
 
        xerr_inside = self.get_value_from_log('X err  inside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_inside, target_val=7.431657E-04, delta=1E-04, label='X err  inside')
 
        xerr_global = self.get_value_from_log('X err   global', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_global, target_val=8.696332E-04, delta=1E-04, label='X err   global')
 
        xerr_outside = self.get_value_from_log('X err  outside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_outside, target_val=8.696332E-04, delta=1E-05, label='X err  outside')
 
        yerr_inside = self.get_value_from_log('Y err  inside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_inside, target_val=9.250194E-04, delta=1E-04, label='Y err  inside')
 
        yerr_global = self.get_value_from_log('Y err   global', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_global, target_val=9.878329E-04, delta=1E-04, label='Y err   global')
 
        yerr_outside = self.get_value_from_log('Y err  outside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_outside, target_val=9.878329E-04, delta=1E-05, label='Y err  outside')
 
        self.assertDelayedFailures()
 
    def tearDown(self):
        self.move_logs()

####################################################################
#  eddy_psi_omega; psi_omega.rea
####################################################################

class Eddy_PsiOmega(NekTestCase):
    example_subdir  = 'eddy_psi_omega'
    case_name        = 'psi_omega'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '16',
            lxd       = '24',
            lx2       = 'lx1-2',
            lelg      = '300',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=1.177007E-10, delta=1E-06, label='X err')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=1.177007E-10, delta=1E-06, label='X err')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

######################################################################
#    rich; eddy_rich.rea
######################################################################

class Eddy_Rich(NekTestCase):
    example_subdir  = 'eddy_rich'
    case_name        = 'eddy_rich'
 
    def setUp(self):
 
        # Default SIZE parameters. Can be overridden in test cases
        self.size_params = dict(
            ldim      = '2',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '500',
        )
 
        self.build_tools(['clean','genmap'])
 
        # Tweak the .rea file and run genmap
        from re import sub
        cls = self.__class__
        rea_path = os.path.join(self.examples_root, cls.example_subdir, cls.case_name + '.rea')
        with open(rea_path, 'r') as f:
            lines = [sub(r'^.*DIVERGENCE$', '      0.10000E-08', l) for l in f]
        with open(rea_path, 'w') as f:
            f.writelines(lines)
        self.run_genmap()
 
    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)
 
        gmres = self.get_value_from_log('gmres ', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=20., label='gmres')
 
        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=5.497982E-05, delta=1E-06, label='X err')
 
        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=8.064398E-05, delta=1E-06, label='Y err')
 
        perr = self.get_value_from_log('P err', column=-5, row=-1)
        self.assertAlmostEqualDelayed(perr, target_val=2.272926E-04, delta=1E-04, label='P err')
 
        self.assertDelayedFailures()
 
    def tearDown(self):
        self.move_logs()

####################################################################
#  eddy; eddy_uv.rea
####################################################################

class Eddy_Uv(NekTestCase):
    example_subdir  = 'eddy_uv'
    case_name        = 'eddy_uv'

    def setUp(self):

        # Default SIZE parameters. Can be overridden in test cases
        self.size_params = dict(
            ldim      = '2',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '500',
        )

        self.build_tools(['clean','genmap'])

        # Tweak the .rea file and run genmap
        from re import sub
        cls = self.__class__
        rea_path = os.path.join(self.examples_root, cls.example_subdir, cls.case_name + '.rea')
        with open(rea_path, 'r') as f:
            lines = [sub(r'^.*DIVERGENCE$', '      0.10000E-08', l) for l in f]
        with open(rea_path, 'w') as f:
            f.writelines(lines)
        self.run_genmap()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-7,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=10., label='gmres')

        crsl = self.get_value_from_log('crsl ', column=-3,)
        self.assertAlmostEqualDelayed(crsl, target_val=0., delta=1500., label='crsl')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=6.007702E-07, delta=1E-08, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=6.489061E-07, delta=1E-08, label='Y err')

        perr = self.get_value_from_log('P err', column=-5, row=-1)
        self.assertAlmostEqualDelayed(perr, target_val=1.448024E-05, delta=1E-06, label='P err')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=8., label='gmres')

        crsl = self.get_value_from_log('crsl ', column=-3,)
        self.assertAlmostEqualDelayed(crsl, target_val=0., delta=1050., label='crsl')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=6.759103E-05, delta=1E-06, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=7.842019E-05, delta=1E-06, label='Y err')

        perr = self.get_value_from_log('P err', column=-5, row=-1)
        self.assertAlmostEqualDelayed(perr, target_val=6.896211E-05, delta=1E-06, label='P err')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  expansion: expansion.rea
####################################################################

class Expansion(NekTestCase):
    example_subdir  = 'expansion'
    case_name        = 'expansion'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '5',
            lxd       = '8',
            lx2       = 'lx1-2',
            lelg      = '2900',
            lpmin     = '2',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
            lfdm      = '0',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek()

        gmres = self.get_value_from_log('gmres', column=-7,row=1)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=80., label='gmres')

        larea = self.get_value_from_log('ubar', column=-4, row=-1)
        self.assertAlmostEqualDelayed(larea, target_val=7.8540E-01, delta=5E-04, label='area')
        
        ubar = self.get_value_from_log('ubar', column=-3, row=-1)
        self.assertAlmostEqualDelayed(ubar, target_val=1.0E-00, delta=5E-04, label='ubar')

        wmax = self.get_value_from_log('ubar', column=-2, row=-1)
        self.assertAlmostEqualDelayed(wmax, target_val=2.0E-00, delta=4E-03, label='ubar')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek()

        gmres = self.get_value_from_log('gmres', column=-6,row=1)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=80., label='gmres')

        larea = self.get_value_from_log('ubar', column=-4, row=-1)
        self.assertAlmostEqualDelayed(larea, target_val=7.8540E-01, delta=5E-04, label='area')
        
        ubar = self.get_value_from_log('ubar', column=-3, row=-1)
        self.assertAlmostEqualDelayed(ubar, target_val=1.0E-00, delta=5E-04, label='ubar')

        wmax = self.get_value_from_log('ubar', column=-2, row=-1)
        self.assertAlmostEqualDelayed(wmax, target_val=2.0E-00, delta=4E-03, label='ubar')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  ext_cyl; ext_cyl.rea
####################################################################

class ExtCyl(NekTestCase):
    example_subdir  = 'ext_cyl'
    case_name        = 'ext_cyl'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '6',
            lxd       = '9',
            lx2       = 'lx1-2',
            lelg      = '1500',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=1000)

        gmres = self.get_value_from_log('gmres', column=-7,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=85., label='gmres')

        dragx = self.get_value_from_log('1dragx', column=-4, row=-1)
        self.assertAlmostEqualDelayed(dragx, target_val=1.2138790E+00, delta=5E-05, label='1dragx')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=1000)

        gmres = self.get_value_from_log('gmres', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=26., label='gmres')

        dragx = self.get_value_from_log('1dragx', column=-4, row=-1)
        self.assertAlmostEqualDelayed(dragx, target_val=1.2139229E+00, delta=5e-05, label='1dragx')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  fs_2; st1.rea, st2.rea, std_wv.rea
####################################################################

class Fs2_St1(NekTestCase):
    example_subdir  = 'fs_2'
    case_name        = 'st1'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '6',
            lxd       = '9',
            lx2       = 'lx1-2',
            lelg      = '100',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=200)

        gmres = self.get_value_from_log('gmres', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=38., label='gmres')

        amp = self.get_value_from_log('amp', column=-2, row=-1)
        self.assertAlmostEqualDelayed(amp, target_val=0.6382379, delta=1e-06, label='amp')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class Fs2_St2(NekTestCase):
    example_subdir  = 'fs_2'
    case_name        = 'st2'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '6',
            lxd       = '9',
            lx2       = 'lx1-2',
            lelg      = '100',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=200)

        gmres = self.get_value_from_log('gmres', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=38., label='gmres')

        amp = self.get_value_from_log('amp', column=-2, row=-1)
        self.assertAlmostEqualDelayed(amp, target_val=0.6376125, delta=1e-06, label='amp')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class Fs2_StdWv(NekTestCase):
    example_subdir  = 'fs_2'
    case_name        = 'std_wv'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '6',
            lxd       = '9',
            lx2       = 'lx1-2',
            lelg      = '100',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=200)

        gmres = self.get_value_from_log('gmres', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=20., label='gmres')

        amp = self.get_value_from_log('amp', column=-2, row=-1)
        self.assertAlmostEqualDelayed(amp, target_val=0.1403036, delta=1e-06, label='amp')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  hemi; hemi
####################################################################

class Hemi(NekTestCase):
    example_subdir = 'hemi'
    case_name = 'hemi'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '5',
            lxd       = '7',
            lx2       = 'lx1',
            lelg      = '2100',
            lpmin     = '2',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-7,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=39., label='gmres')

        wmax = self.get_value_from_log('wmax', column=-2, row=-1)
        self.assertAlmostEqualDelayed(wmax, target_val=4.9173E-01, delta=1e-05, label='wmax')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=34., label='gmres')

        wmax = self.get_value_from_log('wmax', column=-2, row=-1)
        self.assertAlmostEqualDelayed(wmax, target_val=0.48305, delta=1e-05, label='wmax')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  dfh_cav; lin_dfh_cav_dir.par, lin_dfh_cav_adj.par
####################################################################

class LinCav_Adj(NekTestCase):
    example_subdir = 'lin_dfh_cav'
    case_name = 'lin_dfh_cav_adj'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '9',
            lxd       = '13',
            lx2       = 'lx1-2',
            lelg      = '500',
            lpelt     = 'lelt',
        )

        self.build_tools(['clean','genmap'])
        self.run_genmap()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek(usr_file='lin_dfh_cav')
        self.run_nek(step_limit=None)  

        omega = self.get_value_from_log('Energy', column=-3, row=-1)
        self.assertAlmostEqualDelayed(omega, target_val=-7.57304E-03, delta=1E-06, label='growth rate')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class LinCav_Dir(NekTestCase):
    example_subdir = 'lin_dfh_cav'
    case_name = 'lin_dfh_cav_dir'
 
    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '9',
            lxd       = '13',
            lx2       = 'lx1-2',
            lelg      = '500',
            lpelt     = 'lelt',
        )
 
        self.build_tools(['clean','genmap'])
        self.run_genmap()
 
    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek(usr_file='lin_dfh_cav')
        self.run_nek(step_limit=None)  
 
        omega = self.get_value_from_log('Energy', column=-3, row=-1)
        self.assertAlmostEqualDelayed(omega, target_val=-7.57304E-03, delta=1E-06, label='growth rate')
 
        self.assertDelayedFailures()
 
    def tearDown(self):
        self.move_logs()

####################################################################
#  channel2D; lin_chan_dir.par, lin_chan_adj.par
####################################################################

class LinChn_Adj(NekTestCase):
    example_subdir = 'lin_channel2D'
    case_name = 'lin_chan_adj'
 
    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '10',
            lxd       = '15',
            lx2       = 'lx1-2',
            lelg      = '500',
            lpert     = '1',
            lpelt     = 'lelt',
        )
 
        self.build_tools(['clean','genmap'])
        self.run_genmap()
 
    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek(usr_file='lin_chan')
        self.run_nek(step_limit=None)  
 
        omega = self.get_value_from_log('Energy', column=-3, row=-1)
        self.assertAlmostEqualDelayed(omega, target_val=-1.2337E-03, delta=2E-06, label='growth rate')
 
        self.assertDelayedFailures()
 
    def tearDown(self):
        self.move_logs()

class LinChn_Dir(NekTestCase):
    example_subdir = 'lin_channel2D'
    case_name = 'lin_chan_dir'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '10',
            lxd       = '15',
            lx2       = 'lx1-2',
            lelg      = '500',
            lpert     = '1',
            lpelt     = 'lelt',
        )

        self.build_tools(['clean','genmap'])
        self.run_genmap()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek(usr_file='lin_chan')
        self.run_nek(step_limit=None)  

        omega = self.get_value_from_log('Energy', column=-3, row=-1)
        self.assertAlmostEqualDelayed(omega, target_val=-1.2337E-03, delta=2E-06, label='growth rate')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  lowMach_test; lowMach_test.rea
####################################################################

class LowMachTest(NekTestCase):
    example_subdir = 'lowMach_test'
    case_name       = 'lowMach_test'
 
    def setUp(self):
        self.size_params = dict (
            ldim      = '2',
            lx1       = '14',
            lxd       = '20',
            lx2       = 'lx1-0',
            lx1m      = 'lx1',
            lelg      = '500',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap()
 
    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=200)
 
        gmres = self.get_value_from_log(label='gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0, delta=35, label='gmres')
 
        vx = self.get_value_from_log(label='ERROR VX', column=-5, row=-1)
        self.assertAlmostEqualDelayed(vx, target_val=2.6938e-09, delta=1e-10, label='VX')
 
        errt = self.get_value_from_log(label='ERROR T', column=-5, row=-1)
        self.assertAlmostEqualDelayed(errt, target_val=4.5532e-12, delta=1e-13, label='T')
 
        qtl = self.get_value_from_log(label='ERROR QTL', column=-5, row=-1)
        self.assertAlmostEqualDelayed(qtl, target_val=2.6557E-06, delta=1e-07, label='QTL')
 
        self.assertDelayedFailures()
 
    def tearDown(self):
        self.move_logs()


####################################################################
#  mhd; gpf.rea, gpf_m.rea, gpf_b.rea
####################################################################

class Mhd_Gpf(NekTestCase):
    example_subdir = 'mhd'
    case_name = 'gpf'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '6',
            lxd       = '1+3*lx1/2',
            lx2       = 'lx1-2',
            lelg      = '150',
            ldimt     = '2',
            lhis      = '200',
            lelx      = '4',
            lely      = '4',
            lelz      = '8',
            lx1m      = '1',
            lbelt     = 'lelt',
            lpelt     = '1',
            lcvelt    = '1',
            lfdm      = '1',
        )
        self.build_tools(['clean','genbox', 'genmap'])
        self.run_genbox(box_file='gpf')
        self.run_genmap(rea_file='box', tol='0.01')
        self.mvn('box', 'gpf')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        # TODO: This is expected to fail
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        phrase = self.get_phrase_from_log("ABORT: MHD")
        self.assertIsNotNullDelayed(phrase, label='ABORT: MHD')
        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0, delta=15, label='gmres')

        rtavg = self.get_value_from_log('rtavg_gr_Em', column=-4, row=-1)
        self.assertAlmostEqualDelayed(rtavg, target_val=2.56712250E-01, delta=.02, label='rtavg')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()



class Mhd_GpfB(NekTestCase):
    example_subdir = 'mhd'
    case_name = 'gpf_b'

    def setUp(self):
        # Probably a cleaner way to do this...
        # I'm just mimicking the
        self.size_params = dict(
            ldim      = '3',
            lx1       = '6',
            lxd       = '1+3*lx1/2',
            lx2       = 'lx1-2',
            lelg      = '150',
            ldimt     = '2',
            lhis      = '200',
            lelx      = '4',
            lely      = '4',
            lelz      = '8',
            lx1m      = '1',
            lbelt     = 'lelt',
            lpelt     = '1',
            lcvelt    = '1',
            lfdm      = '1',
        )
        self.build_tools(['clean','genbox', 'genmap'])
        self.run_genbox(box_file='gpf')
        self.run_genmap(rea_file='box', tol='0.01')
        self.mvn('box', 'gpf_b')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek(usr_file='gpf')
        self.run_nek(rea_file='gpf_b', step_limit=None)

        phrase = self.get_phrase_from_log("ABORT: MHD")
        self.assertIsNotNullDelayed(phrase, label='ABORT: MHD')
        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek(usr_file='gpf')
        self.run_nek(rea_file='gpf_b', step_limit=None)

        rtavg = self.get_value_from_log('rtavg_gr_Em', column=-4, row=-1)
        self.assertAlmostEqualDelayed(rtavg, target_val=2.56712250E-01, delta=.02, label='rtavg')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class Mhd_GpfM(NekTestCase):
    example_subdir = 'mhd'
    case_name = 'gpf_m'

    def setUp(self):
        import shutil
        # Probably a cleaner way to do this...
        # I'm just mimicking the
        self.size_params = dict(
            ldim      = '3',
            lx1       = '6',
            lxd       = '1+3*lx1/2',
            lx2       = 'lx1-2',
            lelg      = '150',
            ldimt     = '2',
            lhis      = '200',
            lelx      = '4',
            lely      = '4',
            lelz      = '8',
            lx1m      = '1',
            lbelt     = 'lelt',
            lpelt     = '1',
            lcvelt    = '1',
            lfdm      = '1',
        )
        self.build_tools(['clean','genbox', 'genmap'])
        self.run_genbox(box_file='gpf')
        self.run_genmap(rea_file='box', tol='0.01')
        self.mvn('box', 'gpf')
        shutil.copy(
            os.path.join(self.examples_root, self.__class__.example_subdir, 'gpf.map'),
            os.path.join(self.examples_root, self.__class__.example_subdir, 'gpf_m.map')
        )

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek(usr_file='gpf')
        self.run_nek(rea_file='gpf_m', step_limit=None)

        phrase = self.get_phrase_from_log(label="ERROR: FDM")
        self.assertIsNotNullDelayed(phrase, label='ERROR: FDM')
        self.assertDelayedFailures()


    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek(usr_file='gpf')
        self.run_nek(rea_file='gpf_m', step_limit=None)

        rtavg = self.get_value_from_log('rtavg_gr_Em', column=-4, row=-1)
        self.assertAlmostEqualDelayed(rtavg, target_val=2.56712250E-01, delta=.02, label='rtavg')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  moffatt moffat.rea moff_circ.rea
####################################################################

class Moffatt(NekTestCase):
    example_subdir = 'moffatt'
    case_name = 'moffatt'
 
    def setUp(self):
        self.size_params = dict (
            ldim     = '2',
            lx1      = '20',
            lxd      = '22',
            lx2      = 'lx1-2',
            lx1m     = '1',
            lelg     = '100',
            ldimt    = '2',
            lcvelt   = '1',
        )
        self.config_size()
        self.build_tools(['clean','genmap'])
        self.run_genmap()
 
    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek()
 
        phrase = self.get_phrase_from_log('Success: relative strength')
        self.assertIsNotNullDelayed(phrase, label='Success: relative strength')
 
        self.assertDelayedFailures()
 
    def tearDown(self):
        self.move_logs()

class MoffCirc(NekTestCase):
    example_subdir = 'moffatt'
    case_name = 'moff_circ'
 
    def setUp(self):
        self.size_params = dict (
            ldim     = '2',
            lx1      = '20',
            lxd      = '22',
            lx2      = 'lx1-2',
            lx1m     = '1',
            lelg     = '100',
            ldimt    = '2',
            lcvelt   = '1',
        )
        self.config_size()
        self.build_tools(['clean','genmap'])
        self.run_genmap()
 
    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek()
 
        phrase = self.get_phrase_from_log('Success: relative strength')
        self.assertIsNotNullDelayed(phrase, label='Success: relative strength')
 
        self.assertDelayedFailures()
 
    def tearDown(self):
        self.move_logs()

####################################################################
#  mv_cyl with CVODE
####################################################################

class MvCylCvode(NekTestCase):
    example_subdir = 'mv_cyl'
    case_name = 'mv_cyl'
 
    def setUp(self):
        self.size_params = dict (
            ldim     = '2',
            lx1      = '8',
            lxd      = '12',
            lx2      = 'lx1-0',
            lx1m     = 'lx1',
            lelg     = '500',
            ldimt    = '10',
            lcvelt   = 'lelt',
        )
        self.config_size()
        self.build_tools(['clean','genmap'])
        self.run_genmap()
 
    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.log_suffix += '.steps_1e3'
        self.config_parfile({'GENERAL' : {'numSteps' : '1e3', 'dt' : '1e-3'}})
        self.build_nek(opts={'PPLIST':'CVODE'})
        self.run_nek()
 
        err3 = self.get_value_from_log('err', column=-3, row=-1)
        self.assertAlmostEqualDelayed(err3, target_val=0.1743079E-03, delta=1e-6, label='err (column -3)')
 
        err2 = self.get_value_from_log('err', column=-2, row=-1)
        self.assertAlmostEqualDelayed(err2, target_val=0.6348537E-06, delta=1e-9, label='err (column -2)')
 
        self.assertDelayedFailures()
 
    def tearDown(self):
        self.move_logs()

####################################################################
#  mv_wall; mv_wall
####################################################################

class MvWall(NekTestCase):
    example_subdir = 'mv_wall'
    case_name = 'mv_wall'
 
    def setUp(self):
        self.size_params = dict (
            ldim     = '2',
            lx1      = '8',
            lxd      = '12',
            lx2      = 'lx1-0',
            lx1m     = 'lx1',
            lelg     = '500',
            ldimt    = '10',
            lcvelt   = '1',
        )
        self.config_size()
        self.build_tools(['clean','genmap'])
        self.run_genmap()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek()
 
        gmres = self.get_value_from_log(label='gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=47., label='gmres')

        timend = self.get_value_from_log(label='time (non-dimen):', column=-1, row=-1)
        self.assertAlmostEqualDelayed(timend, target_val=1.00, delta=0.00001, label='Final nondimensional time')
 
        self.assertDelayedFailures()
        
    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1 -2'
        self.config_size()
        self.build_nek()
        self.run_nek()
 
        gmres = self.get_value_from_log(label='gmres', column=-6, row=2)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=26., label='gmres')

        timend = self.get_value_from_log(label='time (non-dimen):',column=-1,row=-1)
        self.assertAlmostEqualDelayed(timend, target_val=1.00, delta=0.00001, label='Final nondimensional time')
        
        self.assertDelayedFailures()
 
    def tearDown(self):
        self.move_logs()

####################################################################
#  ocyl: ocyl.rea ocyl2.rea
####################################################################

class Ocyl(NekTestCase):
    example_subdir  = 'ocyl'
    case_name        = 'ocyl'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '10',
            lxd       = '16',
            lx2       = 'lx1-2',
            lelg      = '100',
            lpmin     = '1',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
            lfdm      = '0',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek()

        gmres = self.get_value_from_log('gmres', column=-7,row=1)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=25., label='gmres')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek()

        gmres = self.get_value_from_log('gmres', column=-6,row=1)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=25., label='gmres')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class Ocyl2(NekTestCase):
    example_subdir  = 'ocyl'
    case_name        = 'ocyl2'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '10',
            lxd       = '16',
            lx2       = 'lx1-2',
            lelg      = '100',
            lpmin     = '1',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
            lfdm      = '0',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek()

        gmres = self.get_value_from_log('gmres', column=-7,row=1)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=25., label='gmres')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek()

        gmres = self.get_value_from_log('gmres', column=-6,row=1)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=25., label='gmres')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  os7000; u3_t020_n13.rea
####################################################################

class Os7000(NekTestCase):
    example_subdir = 'os7000'
    case_name = 'u3_t020_n13'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '14',
            lxd       = '3*lx1/2',
            lx2       = 'lx1-2',
            lelg      = '100',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=1000)

        gmres = self.get_value_from_log(label='gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=43., label='gmres')

        egn = self.get_value_from_log(label='egn', column=-2, row=-1)
        self.assertAlmostEqualDelayed(egn, target_val=4.74494769e-05, delta=1e-06, label='egn')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=1000)

        gmres = self.get_value_from_log(label='gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=43., label='gmres')

        egn = self.get_value_from_log(label='egn', column=-2, row=-1)
        self.assertAlmostEqualDelayed(egn, target_val=5.93471252E-05, delta=1e-06, label='egn')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  peris; peris.rea
####################################################################

class Peris(NekTestCase):
    example_subdir = 'peris'
    case_name = 'peris'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '6',
            lxd       = '9',
            lx2       = 'lx1',
            lelg      = '200',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log(label='gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=26., label='gmres')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log(label='gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=18., label='gmres')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  pipe; helix.rea, stenosis.rea
####################################################################

class Pipe_Helix(NekTestCase):
    example_subdir = 'pipe'
    case_name = 'helix'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '6',
            lxd       = '10',
            lx2       = 'lx1-2',
            lelg      = '500',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=61., label='gmres')

        err2 = self.get_value_from_log('err2', column=-2, row=-1)
        self.assertAlmostEqualDelayed(err2, target_val=1.9077617E+00, delta=1e-06, label='err2')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=123., label='gmres')

        err2 = self.get_value_from_log('err2', column=-2, row=-1)
        self.assertAlmostEqualDelayed(err2, target_val=1.9072258E+00, delta=1e-06, label='err2')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class Pipe_Stenosis(NekTestCase):
    example_subdir = 'pipe'
    case_name = 'stenosis'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '6',
            lxd       = '10',
            lx2       = 'lx1-2',
            lelg      = '500',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        n2to3_input = [
            'w2dcyl020a',
            'stenosis',
            '0                      ascii output',
            '20                     input number of levels: (1, 2, 3,... etc.?):',
            '0                      input z min:',
            '10                     input z max:',
            '1                      input gain (0=custom,1=uniform,other=geometric spacing):',
            'n                      This is for CEM: yes or no:',
            'v                      Enter Z (5) boundary condition (P,v,O):',
            'O                      Enter Z (6) boundary condition (v,O):',
            'y                      Formatted .rea file? (y or Y):',
        ]
        self.build_tools(['clean','n2to3','genmap'])
        self.run_n2to3(n2to3_input)
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=196., label='gmres')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=51., label='gmres')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  rayleigh; ray1.rea, ray2.rea
####################################################################

class Rayleigh_Ray1(NekTestCase):
    example_subdir  = 'rayleigh'
    case_name        = 'ray1'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '12',
            lxd       = '18',
            lx2       = 'lx1',
            lelg      = '50',
            ldimt     = '1',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(rea_file='ray1', tol='0.01')


    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek(usr_file='ray0')
        self.run_nek(rea_file='ray1', step_limit=200)

        gmres = self.get_value_from_log(label='gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=32., label='gmres')

        umax = self.get_value_from_log(label='umax', column=-3, row=-1)
        self.assertAlmostEqualDelayed(umax, target_val=0.00377251, delta=2e-05, label='umax')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek(usr_file='ray0')
        self.run_nek(rea_file='ray1', step_limit=200)

        gmres = self.get_value_from_log(label='gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=11., label='gmres')

        umax = self.get_value_from_log(label='umax', column=-3, row=-1)
        self.assertAlmostEqualDelayed(umax, target_val=4.831113E-03, delta=2e-05, label='umax')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class Rayleigh_Ray2(NekTestCase):
    example_subdir  = 'rayleigh'
    case_name        = 'ray2'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '12',
            lxd       = '18',
            lx2       = 'lx1',
            lelg      = '50',
            ldimt     = '1',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genbox','genmap'])
        self.run_genbox(box_file='ray2')
        self.run_genmap(rea_file='box', tol='0.01')
        self.mvn('box', 'ray2')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek(usr_file='ray0')
        self.run_nek(rea_file='ray2', step_limit=200)

        gmres = self.get_value_from_log('gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0, delta=31, label='gmres')

        umax = self.get_value_from_log('umax', column=-3, row=-1)
        self.assertAlmostEqualDelayed(umax, target_val=0.00551992, delta=1e-04, label='umax')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek(usr_file='ray0')
        self.run_nek(rea_file='ray2', step_limit=200)

        gmres = self.get_value_from_log('gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0, delta=11, label='gmres')

        umax = self.get_value_from_log('umax', column=-3, row=-1)
        self.assertAlmostEqualDelayed(umax, target_val=0.006621465, delta=1e-04, label='umax')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  expansion: expansion.rea
####################################################################

class Robin(NekTestCase):
    example_subdir  = 'robin'
    case_name        = 'robin'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '12',
            lxd       = '18',
            lx2       = 'lx1-0',
            lelg      = '10',
            lpmin     = '1',
            lpmax     = '1',
            ldimt     = '3',
            lhis      = '100',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
            lfdm      = '0',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_serial
    def test_PnPn_Serial(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek()

        hmh = self.get_value_from_log('Hmholtz TEMP', column=-4,row=-1)
        self.assertAlmostEqualDelayed(hmh, target_val=0., delta=80., label='Hmholtz TEMP')

        tbar = self.get_value_from_log('tbar', column=-2, row=-1)
        self.assertAlmostEqualDelayed(tbar, target_val=0.0E-00, delta=1E-05, label='tbar')

        self.assertDelayedFailures()

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek()

        hmh = self.get_value_from_log('Hmholtz TEMP', column=-4,row=-1)
        self.assertAlmostEqualDelayed(hmh, target_val=0., delta=20., label='Hmholtz TEMP')

        tbar = self.get_value_from_log('tbar', column=-2, row=-1)
        self.assertAlmostEqualDelayed(tbar, target_val=0.0E-00, delta=1E-05, label='tbar')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  shear4; shear4.rea, thin.rea
####################################################################

class Shear4_Shear4(NekTestCase):
    example_subdir = 'shear4'
    case_name = 'shear4'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '9',
            lxd       = '13',
            lx2       = 'lx1-2',
            lelg      = '300',
            ldimt     = '4',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '16',
            lely      = '16',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=26., label='gmres')

        vort = self.get_value_from_log('peak vorticity', column=-3, row=-1)
        self.assertAlmostEqualDelayed(vort, target_val=3.031328E+01, delta=1e-06, label='peak vorticity')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=17., label='gmres')

        vort = self.get_value_from_log('peak vorticity', column=-3, row=-1)
        self.assertAlmostEqualDelayed(vort, target_val=3.031328E+01, delta=1e-06, label='peak vorticity')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class Shear4_Thin(NekTestCase):
    example_subdir = 'shear4'
    case_name = 'thin'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '9',
            lxd       = '13',
            lx2       = 'lx1-2',
            lelg      = '300',
            ldimt     = '4',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '16',
            lely      = '16',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=26., label='gmres')

        vort = self.get_value_from_log('peak vorticity', column=-3, row=-1)
        self.assertAlmostEqualDelayed(vort, target_val=9.991753E+01, delta=1e-06, label='peak vorticity')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=17., label='gmres')

        vort = self.get_value_from_log('peak vorticity', column=-3, row=-1)
        self.assertAlmostEqualDelayed(vort, target_val=9.991556E+01, delta=1e-06, label='peak vorticity')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()


####################################################################
#  smoother; lpt.par
####################################################################

class Smooth(NekTestCase):
    example_subdir = 'smoother'
    case_name = 'lpt'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '3',
            lxd       = '4',
            lx2       = 'lx1',
            lelg      = '10000',
            ldimt     = '1',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        phrase = self.get_phrase_from_log('Mesh smoothing completed.')
        self.assertIsNotNullDelayed(phrase, label='Mesh smoothing completed.')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  solid; solid.rea
####################################################################

class Solid(NekTestCase):
    example_subdir = 'solid'
    case_name = 'solid'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '50',
            ldimt     = '1',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        error = self.get_value_from_log('error', column=-2, row=-1)
        self.assertAlmostEqualDelayed(error, target_val=7.821228E-05, delta=1e-06, label='error')
        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        error = self.get_value_from_log('error', column=-2, row=-1)
        self.assertAlmostEqualDelayed(error, target_val=7.821228E-05, delta=1e-06, label='error')
        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  strat; re10f1000p1000.rea, re10f1000p0001.rea
####################################################################

class Strat_P0001(NekTestCase):
    example_subdir = 'strat'
    case_name = 're10f1000p0001'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '400',
            ldimt     = '1',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=200)

        gmres = self.get_value_from_log('gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0, delta=60, label='gmres')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=200)

        upres = self.get_value_from_log('U-PRES', column=-6)
        self.assertAlmostEqualDelayed(upres, target_val=0, delta=27, label='U-PRES')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class Strat_P1000(NekTestCase):
    example_subdir = 'strat'
    case_name = 're10f1000p1000'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '400',
            ldimt     = '1',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=200)

        gmres = self.get_value_from_log('gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0, delta=60, label='gmres')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=200)

        upres = self.get_value_from_log('U-PRES', column=-6)
        self.assertAlmostEqualDelayed(upres, target_val=0, delta=27, label='U-PRES')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  taylor; taylor.rea
####################################################################

class Taylor(NekTestCase):
    example_subdir = 'taylor'
    case_name = 'taylor'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '10',
            lxd       = '16',
            lx2       = 'lx1-2',
            lelg      = '20',
            ldimt     = '1',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=23., label='gmres')

        tq = self.get_value_from_log('tq', column=-5, row=-1)
        self.assertAlmostEqualDelayed(tq, target_val=4.13037E-06, delta=1e-06, label='tq')

        err = self.get_value_from_log('err', column=-2, row=-1)
        self.assertAlmostEqualDelayed(err, target_val=2.973648E-09, delta=1e-06, label='err')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        solver_time = self.get_value_from_log('total solver time', column=-2)
        self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=40., label='total solver time')

        gmres = self.get_value_from_log('gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=14, label='gmres')

        tq = self.get_value_from_log('tq', column=-5, row=-1)
        self.assertAlmostEqualDelayed(tq, target_val=4.10783E-06, delta=1e-06, label='tq')

        err = self.get_value_from_log('err', column=-2, row=-1)
        self.assertAlmostEqualDelayed(err, target_val=2.826284E-10, delta=1e-06, label='err')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

###############################################################################
#  turbChannel: turbChannel.rea
###############################################################################

class TurbChannel(NekTestCase):
    example_subdir = 'turbChannel'
    case_name = 'turbChannel'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '1600',
            lpmin     = '4',
            ldimt     = '3',
            lhis      = '100',
            lelx      = '8',
            lely      = '8',
            lelz      = '8',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=35., label='gmres')

#        solver_time = self.get_value_from_log('total solver time', column=-2)
#        self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=37.0, label='total solver time')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=35., label='gmres')

#        solver_time = self.get_value_from_log('total solver time', column=-2)
#        self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=25.0, label='total solver time')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

###############################################################################
#  turbInflow: turbInflow.rea
###############################################################################

class TurbInflow(NekTestCase):
    example_subdir = 'turbInflow'
    case_name = 'turbInflow'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '700',
            lpmin     = '1',
            ldimt     = '1',
            lhis      = '1000',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        from subprocess import call, check_call, Popen, PIPE, STDOUT
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=40)

        gmres = self.get_value_from_log('gmres', column=-7)
        e2 = self.get_value_from_log('e2', column=-4, row=-1)
        ub = self.get_value_from_log('e2', column=-3, row=-1)

        cwd = os.path.join(self.examples_root, self.example_subdir)
        genlist_in  = os.path.join(cwd, 'genlist')
        proc = Popen([genlist_in],cwd=cwd,shell=True,stdin=PIPE)
        proc.wait()
        
        # move log file to avoid it to be ovewritten
        for f in os.listdir(cwd):
                if 'log' in f and '_run' not in f:
                    src_file = os.path.join(cwd, f)
                    dest_file = os.path.join(cwd, f+'_run')
                    try:
                        os.rename(src_file, dest_file)
                    except OSError as E:
                        # TODO: change to warnings.warning
                        print("    Could not move {0} to {1}: {2}".format(src_file, dest_file, E))
                    else:
                        print("    Moved {0} to {1}".format(src_file, dest_file))
        
        self.build_nek(usr_file='prms')
        self.run_nek(step_limit=None)
        phrase = self.get_phrase_from_log("FILE:avxturbInflow0.f00001")

        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=45., label='gmres')
        self.assertAlmostEqualDelayed(e2, target_val=7.6E-03, delta=2.0E-04, label='e2')
        self.assertAlmostEqualDelayed(ub, target_val=1.0250, delta=2.0E-04, label='ub')
        self.assertIsNotNullDelayed(phrase, label='FILE:avxturbInflow0.f00001')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        from subprocess import call, check_call, Popen, PIPE, STDOUT
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=40)

        gmres = self.get_value_from_log('gmres', column=-6)
        e2 = self.get_value_from_log('e2', column=-4, row=-1)
        ub = self.get_value_from_log('e2', column=-3, row=-1)

        cwd = os.path.join(self.examples_root, self.example_subdir)
        genlist_in  = os.path.join(cwd, 'genlist')
        proc = Popen([genlist_in],cwd=cwd,shell=True,stdin=PIPE)
        proc.wait()
        
        # move log file to avoid it to be ovewritten
        for f in os.listdir(cwd):
                if 'log' in f and '_run' not in f:
                    src_file = os.path.join(cwd, f)
                    dest_file = os.path.join(cwd, f+'_run')
                    try:
                        os.rename(src_file, dest_file)
                    except OSError as E:
                        # TODO: change to warnings.warning
                        print("    Could not move {0} to {1}: {2}".format(src_file, dest_file, E))
                    else:
                        print("    Moved {0} to {1}".format(src_file, dest_file))
        
        self.build_nek(usr_file='prms')
        self.run_nek(step_limit=None)
        phrase = self.get_phrase_from_log("FILE:avxturbInflow0.f00001")

        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=90., label='gmres')
        self.assertAlmostEqualDelayed(e2, target_val=7.6E-03, delta=2.0E-04, label='e2')
        self.assertAlmostEqualDelayed(ub, target_val=1.0250, delta=2.0E-04, label='ub')
        self.assertIsNotNullDelayed(phrase, label='FILE:avxturbInflow0.f00001')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  vortex; r1854a.rea
####################################################################

class Vortex(NekTestCase):
    example_subdir = 'vortex'
    case_name = 'r1854a'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '200',
            ldimt     = '2',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=65., label='gmres')

        vmin = self.get_value_from_log('VMIN', column=-2, row=-1)
        self.assertAlmostEqualDelayed(vmin, target_val=-1.910312E-03, delta=1e-05, label='VMIN')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=10)

        gmres = self.get_value_from_log('gmres', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=18., label='gmres')

        vmin = self.get_value_from_log('VMIN', column=-2, row=-1)
        self.assertAlmostEqualDelayed(vmin, target_val=-1.839120E-03, delta=1e-05, label='VMIN')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  vortex2; v2d
####################################################################

class Vortex2(NekTestCase):
    example_subdir = 'vortex2'
    case_name = 'v2d'

    def setUp(self):
        import re
        self.size_params = dict(
            ldim      = '2',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '300',
            ldimt     = '4',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '20',
            lely      = '60',
            lelz      = '1',
            mprev     = '80',
            lgmres    = '40',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['clean','genmap'])
        self.run_genmap(tol='0.01')

        # Tweak .rea file
        rea_file_path = os.path.join(self.examples_root, self.__class__.example_subdir, self.case_name + '.rea')
        with open(rea_file_path, 'r') as f:
            lines = [re.sub(r'(^\s+[\d.]+\s+p11.*$)', r' 8000\g<1>', l) for l in f]
        with open(rea_file_path, 'w') as f:
            f.writelines(lines)

        # Extra tweaks to the SIZE file
        size_file_path = os.path.join(self.examples_root, self.__class__.example_subdir, 'SIZE')
        with open(size_file_path, 'r') as f:
            lines = f.readlines()

        lines = [re.sub(
            r'( {6}parameter *)\(lx1=10,ly1=lx1,lz1=1,lelt=80,lelv=lelt\)( *)',
            r'\g<1>(lx1=8,ly1=lx1,lz1=1,lelt=80,lelv=lelt)\g<2>', l, flags=re.I) for l in lines]
        lines = [re.sub(
            r'( {6}parameter *)\(lxd=15,lyd=lxd,lzd=1\)( *)',
            r'\g<1>(lxd=12,lyd=lxd,lzd=1)\g<2>', l, flags=re.I) for l in lines]

        with open(size_file_path, 'w') as f:
            f.writelines(lines)

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        pres = self.get_value_from_log('PRES', column=-4)
        self.assertAlmostEqualDelayed(pres, target_val=0., delta=100., label='PRES')

        umin = self.get_value_from_log('umin', column=-2, row=-1)
        self.assertAlmostEqualDelayed(umin, target_val=-1.453402E-03, delta=1e-03, label='umin')

        torqx = self.get_value_from_log('1torqx', column=-2, row=-1)
        self.assertAlmostEqualDelayed(torqx, target_val=-1.7399905E-07, delta=1e-06, label='1torqx')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        upress = self.get_value_from_log('U-PRES', column=-5)
        self.assertAlmostEqualDelayed(upress, target_val=0., delta=100, label='U-PRES')

        umin = self.get_value_from_log('umin', column=-2, row=-1)
        self.assertAlmostEqualDelayed(umin, target_val=-2.448980E-03, delta=1e-03, label='umin')

        torqx = self.get_value_from_log('1torqx', column=-2, row=-1)
        self.assertAlmostEqualDelayed(torqx, target_val=-1.6276138E-07, delta=1e-06, label='1torqx')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()
        
###############################################################

if __name__ == '__main__':
    import unittest, argparse, os

    # Get arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--f77", default='mpif77', help="The Fortran 77 compiler to use [default: mpif77]")
    parser.add_argument("--cc", default='mpicc',  help="The C compiler to use [default: mpicc]")
    parser.add_argument("--ifmpi", default='true', choices=['true', 'false'], help="Enable/disable parallel tests with MPI [default: true]")
    parser.add_argument("--nprocs", default='4', help="Number of processes to use for MPI tests [default: 4]")
    parser.add_argument("-v", "--verbose", action='store_true', help="Enable verbose output")
 
    args = parser.parse_args()

    # # Set environment
    os.environ['CC'] = args.cc
    os.environ['FC'] = args.f77
    os.environ['IFMPI'] = args.ifmpi
    os.environ['PARALLEL_PROCS'] = args.nprocs
    if args.verbose:
        os.environ['VERBOSE_TESTS'] = 'true'
        ut_verbose = 2
    else:
        os.environ['VERBOSE_TESTS'] = 'false'
        ut_verbose = 1

    testList = (
               Annulus2d,
               Benard_Ray9,
               Blasius,
               CmtInviscidVortex,
               Cone16,
               Cone64,
               Cone256,
               ConjHt,
               CylRestart_C,
               CylRestart_P,
               Eddy_Neknek,
               Eddy_PsiOmega, 
               Eddy_Rich,
               Eddy_Uv,
               Expansion,
               ExtCyl,
               Fs2_St1,
               Fs2_St2,
               Fs2_StdWv,
               Hemi,
               LinCav_Adj,
               LinCav_Dir,
               LinChn_Adj,
               LinChn_Dir,
               LowMachTest,
               Mhd_Gpf,
               Mhd_GpfB,
               Mhd_GpfM,
               Moffatt,
               MoffCirc,
               MvCylCvode,
               MvWall,
               Ocyl,
               Ocyl2,
               Os7000,
               Peris,
               Pipe_Helix,
               Pipe_Stenosis,
               Rayleigh_Ray1,
               Rayleigh_Ray2,
               Robin,
               Shear4_Shear4,
               Shear4_Thin,
               Smooth,
               Solid,
               Strat_P0001,
               Strat_P1000,
               Taylor,
               TurbChannel,
               TurbInflow,
               Vortex,
               Vortex2
               ) 

    suite = unittest.TestSuite([unittest.TestLoader().loadTestsFromTestCase(t) for t in testList])
    unittest.TextTestRunner(verbosity=ut_verbose, buffer=True).run(suite)

