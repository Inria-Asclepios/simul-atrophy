#!/usr/bin/env python
import os
import os.path as op
import sys
import optparse
import bish_utils as bu

def get_input_options():
    """ command line interface, get input options and interact with the user
    """
    #parse the command line inputs
    parser = optparse.OptionParser()
    parser.add_option('-p', '--patientID', dest='patient', help='patient dir at'
                      ' $ADLEM_DIR/results/patients/ which contains all the '
                      'inputs required by this script. Results will also be '
                      ' written here.')
    parser.add_option('-a', '--atrophy', dest='atrophy', help='input atrophy '
                      'file.')
    parser.add_option('-s', '--in_seg', dest='in_seg', help='input brain seg ')
    parser.add_option('--use_dti', dest='use_dti', action='store_true', help=
                      'If the option is provided, input DTI must be given with'
                      ' the option -d')
    parser.add_option('-d', '--in_dti', dest='in_dti', help='valid '
                      'tensor image file present in patientID directory.')
    parser.add_option('-f', '--res_prefix', dest='res_prefix', help='prefix to'
                      ' be added to all the output filenames')
    parser.add_option('-i', '--in_img', dest='in_img', help='valid input image '
                      'that will be warped by the obtained displacement fields '
                      'from the model.')
    parser.add_option('-t', '--time_steps', dest='time_steps', help='(integer) '
                      'number of time steps you want to solve the system.')
    parser.add_option('--wrt_press', dest='wrt_press', action='store_true',
                      help='Write pressure map on the disk.')
    parser.add_option('--wrt_force', dest='wrt_force', action='store_true',
                      help='Write force on the the disk.')
    parser.add_option('--wrt_residual', dest='wrt_residual', action='store_true'
                      , help='Write residual on the disk.')
    parser.add_option('--in_cluster', dest='in_cluster', action='store_true',
                      help='If the option is provided, launches as a job on the'
                      ' cluster.')

    (ops, args) = parser.parse_args()
    if ops.patient is None:
        ops.patient = raw_input('Enter a valid existing patient directory name')
    if ops.time_steps is None:
        ops.time_steps = raw_input('Enter num_of_time_steps to solve.')
    if ops.in_img is None:
        ops.in_img = raw_input('Enter image file that will be warped by the '
                               'obtained displacement fields.')
    if ops.atrophy is None:
        ops.atrophy = raw_input('Enter a valid atrophy file.')
    if ops.in_seg is None:
        ops.in_seg = raw_input('Enter a valid segmented image ')
    if ops.use_dti is True:
        if ops.in_dti is None:
            ops.in_dti = raw_input('Enter a DTI file (since you are using '
                                   '--use_dti.):  ')
    if ops.res_prefix is None:
        ops.res_prefix = raw_input('Enter a prefix string to be added to all '
                                   'the output filenames')
    return ops


def get_petsc_options():
    """ Return petsc options. Change the code as required for different
    options you want. Always puts one character space at the beginning and at
    the end.
    """
    # Put one character space at the beginning and the end!
    ops = (' -pc_type fieldsplit -pc_fieldsplit_type schur '
           '-pc_fieldsplit_schur_precondition self -pc_fieldsplit_dm_splits 0 '
           '-pc_fieldsplit_0_fields 0,1,2 -pc_fieldsplit_1_fields 3 '
           '-fieldsplit_0_pc_type hypre ')
    #monitor options
    # ops = ops + ('-fieldsplit_0_ksp_converged_reason -fieldsplit_0_ksp_max_it'
    #            ' 100 -ksp_converged_reason')
    ops = ops + ('-fieldsplit_1_ksp_converged_reason -ksp_converged_reason ')
    # ops = ops + ('-fieldsplit_1_ksp_max_it 3 -ksp_max_it 3 -ksp_rtol 1.0e-8 '
    #              '-ksp_monitor_true_residual')
    ops = ops + ('-fieldsplit_1_ksp_monitor_true_residual -ksp_monitor_true_'
                 'residual -log_summary -ksp_view ')
    return ops


def main():
    """
    Run the model adlem for a given patient, segmentation mask, atrophy map and
    [optional DTI]."""

    work_dir = os.getenv('ADLEM_DIR')
    if work_dir is None:
        raise ValueError('environment variable ADLEM_DIR not set.')
    target = op.join(work_dir, 'build/src/AdLemMain')
    ops = get_input_options()
    res_dir = op.join(work_dir, 'results/patients', ops.patient)
    in_img = op.join(res_dir, ops.in_img)
    atrophy = op.join(res_dir, ops.atrophy)
    in_seg = op.join(res_dir, ops.in_seg)
    res_path = res_dir+'/' #due to stupid limitation in my AdLem code
    # where it concatenates simply the -resPath argument with the file prefix.
    if ops.use_dti is True:
        use_dti = 'true'
        in_dti = op.join(res_dir, ops.in_dti)
    else:
        use_dti = 'false'
        in_dti = "dummyFileName"
    if ops.wrt_press is True:
        wrt_press = 'true'
    else:
        wrt_press = 'false'
    if ops.wrt_force is True:
        wrt_force = 'true'
    else:
        wrt_force = 'false'
    if ops.wrt_residual is True:
        wrt_residual = 'true'
    else:
        wrt_residual = 'false'
    petsc_ops = get_petsc_options()
    cmd = target + ' -atrophyFile ' + atrophy + ' -useTensorLambda ' + use_dti+\
          ' -lambdaFile ' + in_dti + ' -maskFile ' + in_seg + ' -imageFile ' + \
          in_img + ' -numOfTimeSteps ' + ops.time_steps + ' -resPath '+res_path\
          + ' -resultsFilenamesPrefix ' + ops.res_prefix + ' -writePressure ' +\
          wrt_press + ' -writeForce ' + wrt_force + ' -writeResidual ' + \
          wrt_residual + petsc_ops
    if ops.in_cluster is True:
        cluster_mpi = '/opt/openmpi-gcc/current/bin/mpiexec '
        cmd = bu.sophia_nef_pbs_setting() + cluster_mpi + cmd
        job_name = 'AdLemModel' + ops.patient
        #print cmd
        bu.qsub_job(name=job_name, queue='parlong', walltime='12:00:00',
                    procs='nodes=3:xeon:ppn=20', mem='mem=180gb',
                    dest_dir=res_dir, cmd=cmd)
    else:
        #print cmd
        bu.print_and_execute(cmd)


if __name__ == "__main__":
    main()
