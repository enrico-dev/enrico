from lib.nekTestCase import *
from lib.nekBinBuild import build_nek
from lib.nekBinRun import run_nek
import re

class AccTestCase(NekTestCase):

    step_limit = None
    log_label = ''

    @staticmethod
    def get_lines_from_log(label, logfile):
        with open(logfile, 'r') as f:
            line_list = [l for l in f if re.search(r'\b{0}\b'.format(label), l)]
        return line_list

    def compare_runs(self):

        cls = self.__class__

        # First get the output w/o acc

        build_nek(
            source_root = self.source_root,
            usr_file    = cls.case_name,
            cwd         = os.path.join(self.examples_root, cls.example_subdir),
            f77         = self.f77.replace('-acc',''),
            cc          = self.cc.replace('-acc',''),
            ifmpi       = self.ifmpi,
            verbose     = self.verbose
        )

        run_nek(
            cwd        = os.path.join(self.examples_root, cls.example_subdir),
            rea_file   = cls.case_name,
            ifmpi      = self.ifmpi,
            log_suffix = self.log_suffix + '.noacc',
            n_procs    = self.mpi_procs,
            step_limit = cls.step_limit,
            verbose    = self.verbose
        )

        # Then get the output with  ACC

        build_nek(
            source_root = self.source_root,
            usr_file    = cls.case_name,
            cwd         = os.path.join(self.examples_root, cls.example_subdir),
            f77         = self.f77.replace('-acc','') + ' -acc -Minfo=accel',
            cc          = self.cc.replace('-acc','') + ' -acc -Minfo=accel',
            ifmpi       = self.ifmpi,
            verbose     = self.verbose
        )

        run_nek(
            cwd        = os.path.join(self.examples_root, cls.example_subdir),
            rea_file   = cls.case_name,
            ifmpi      = self.ifmpi,
            log_suffix = self.log_suffix + '.acc',
            n_procs    = self.mpi_procs,
            step_limit = cls.step_limit,
            verbose    = self.verbose
        )

        # Now compare output

        logfile_noacc = os.path.join(
            self.examples_root,
            cls.example_subdir,
            '{0}.log.{1}{2}.noacc'.format(cls.case_name, self.mpi_procs, self.log_suffix)
        )

        logfile_acc = os.path.join(
            self.examples_root,
            cls.example_subdir,
            '{0}.log.{1}{2}.acc'.format(cls.case_name, self.mpi_procs, self.log_suffix)
        )

        loglines_noacc = self.get_lines_from_log(cls.log_label, logfile_noacc)
        loglines_acc = self.get_lines_from_log(cls.log_label, logfile_acc)

        self.assertNotEqual(
            len(logfile_noacc),
            0,
            'Could not find the label {0} in the logfile {1}'.format(cls.log_label, logfile_noacc)
        )

        self.assertNotEqual(
            len(loglines_acc),
            0,
            'Could not find the label {0} in the logfile {1}'.format(cls.log_label, logfile_acc)
        )

        for (noacc, acc) in zip(loglines_noacc, loglines_acc):
            err_msg = 'The lines with/without ACC did not match after this point\n' + \
                      '\n---- Without ACC -----------------------------------------------------\n' + \
                      noacc + \
                      '\n---- With ACC --------------------------------------------------------\n' + \
                      acc
            self.assertEqual(noacc, acc, err_msg)

#===================================================================
#  pipe; stenosis.rea
#===================================================================

class Compare_Pipe_Stenosis(AccTestCase):
    example_subdir = 'pipe'
    case_name = 'stenosis'
    step_limit = 10
    log_label = 'DIV'

    def setUp(self):
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
        self.build_tools(['n2to3', 'genmap'])
        self.run_n2to3(n2to3_input)
        self.run_genmap()

    @pn_pn_serial
    def test_PnPn_Serial(self):
        self.config_size(lx2='lx1', ly2='ly1', lz2='lz1')
        self.compare_runs()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.config_size(lx2='lx1', ly2='ly1', lz2='lz1')
        self.compare_runs()

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
        self.config_size(lx2='lx1-2', ly2='ly1-2', lz2='lz1-2')
        self.compare_runs()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.config_size(lx2='lx1-2', ly2='ly1-2', lz2='lz1-2')
        self.compare_runs()

    def tearDown(self):
        self.move_logs()

#===================================================================
#  solid; solid.rea
#===================================================================

class Compare_Solid(AccTestCase):
    example_subdir = 'solid'
    case_name = 'solid'
    step_limit = None
    log_label = 'error'

    def setUp(self):
        self.build_tools(['genmap'])
        self.run_genmap()

    @pn_pn_serial
    def test_PnPn_Serial(self):
        self.config_size(lx2='lx1', ly2='ly1', lz2='lz1')
        self.compare_runs()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.config_size(lx2='lx1', ly2='ly1', lz2='lz1')
        self.compare_runs()

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
        self.config_size(lx2='lx1-2', ly2='ly1-2', lz2='lz1-2')
        self.compare_runs()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.config_size(lx2='lx1-2', ly2='ly1-2', lz2='lz1-2')
        self.compare_runs()

    def tearDown(self):
        self.move_logs()

#===================================================================
#  3dbox: b3d.rea
#===================================================================

class Compare_ThreeDBox(AccTestCase):
    example_subdir = '3dbox'
    case_name = 'b3d'
    step_limit = 10
    log_label = 'DIV'

    def setUp(self):
        self.build_tools(['genbox', 'genmap'])
        self.run_genbox()
        self.mvn('box', self.__class__.case_name)
        self.run_genmap()

    @pn_pn_serial
    def test_PnPn_Serial(self):
        self.config_size(lx2='lx1', ly2='ly1', lz2='lz1')
        self.compare_runs()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.config_size(lx2='lx1', ly2='ly1', lz2='lz1')
        self.compare_runs()

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
        self.config_size(lx2='lx1-2', ly2='ly1-2', lz2='lz1-2')
        self.compare_runs()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.config_size(lx2='lx1-2', ly2='ly1-2', lz2='lz1-2')
        self.compare_runs()

    def tearDown(self):
        self.move_logs()

#===================================================================
#  turbChannel; turbChannel.rea
#===================================================================

class Compare_TurbChannel(AccTestCase):
    example_subdir = 'turbChannel'
    case_name = 'turbChannel'
    step_limit = 10
    log_label = 'err2'

    def setUp(self):
        self.build_tools(['genmap'])
        self.run_genmap(tol='0.5')

    @pn_pn_serial
    def test_PnPn_Serial(self):
        self.config_size(lx2='lx1', ly2='ly1', lz2='lz1')
        self.compare_runs()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.config_size(lx2='lx1', ly2='ly1', lz2='lz1')
        self.compare_runs()

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
        self.config_size(lx2='lx1-2', ly2='ly1-2', lz2='lz1-2')
        self.compare_runs()

    @pn_pn_2_parallel
    def test_PnPn2_Serial(self):
        self.config_size(lx2='lx1-2', ly2='ly1-2', lz2='lz1-2')
        self.compare_runs()

    def tearDown(self):
        self.move_logs()

#===================================================================
#  vortex; r1854a.rea
#===================================================================

class Compare_Vortex(AccTestCase):
    example_subdir = 'vortex'
    case_name = 'r1854a'
    step_limit = 10
    log_label = 'VMIN'

    def setUp(self):
        self.build_tools(['genmap'])
        self.run_genmap()

    @pn_pn_serial
    def test_PnPn_Serial(self):
        self.config_size(lx2='lx1', ly2='ly1', lz2='lz1')
        self.compare_runs()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.config_size(lx2='lx1', ly2='ly1', lz2='lz1')
        self.compare_runs()

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
        self.config_size(lx2='lx1-2', ly2='ly1-2', lz2='lz1-2')
        self.compare_runs()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.config_size(lx2='lx1-2', ly2='ly1-2', lz2='lz1-2')
        self.compare_runs()

    def tearDown(self):
        self.move_logs()
