#!/usr/bin/env python
import os
import os.path as op
import argparse as ag
import bish_utils as bu

def get_input_options():
    """ command line interface, get input options and interact with the user
    """
    #parse the command line inputs
    parser = ag.ArgumentParser()
    parser.add_argument(
        'patient', help='patient dir at $ADLEM_DIR/results/patients/ which '
        'contains all the inputs required by this script. Results will also be '
        ' written here.')
    parser.add_argument(
        'res_prefix', help='prefix to be added to all the output filenames')
    parser.add_argument(
        'lame_paras', help='input lame parameters: Format: muTissue,muCsf,'
        'lambdaTissue,lambdaCsf \n If both --mu_image and --in_dti used, '
        'lame_paras wont be used and those images will be used instead.')
    parser.add_argument('boundary_condition', help='Possible values: '
                        'dirichlet_at_walls dirichlet_at_skull')
    parser.add_argument('atrophy', help='input atrophy file.')
    parser.add_argument('in_seg', help='input brain seg ')
    parser.add_argument(
        'in_img', help='valid input image that will be warped by the obtained'
        ' displacement fields from the model.')
    parser.add_argument('petsc_op_file', help='file with petsc options')
    parser.add_argument('time_steps', help='(integer) '
                        'number of time steps you want to solve the system.')
    parser.add_argument('--div12pt_stencil', action='store_true',
                        help='If given, use 12pt stencil for divergence.')
    parser.add_argument(
        '--relax_ic_in_csf', action='store_true', help='If given, relaxes IC.')
    parser.add_argument(
        '-k', '--relax_ic_coeff', help='relaxation coeff. Relevant only when '
        'using --relax_ic_in_csf. If --relax_ic_in_csf given and this option is'
        ' not provided, uses reciprocal of lamda as coefficient.')
    parser.add_argument('-m', '--mu_file', help='input mu image file: ')
    parser.add_argument(
        '--use_dti', dest='use_dti', action='store_true', help='If the option '
        'is provided, input DTI must be given with the option -d')
    parser.add_argument(
        '-d', '--in_dti', help='tensor image file in patientID directory.')
    parser.add_argument(
        '--wrt_press', action='store_true', help='Write pressure file.')
    parser.add_argument(
        '--wrt_force', action='store_true', help='Write force file.')
    parser.add_argument(
        '--wrt_residual', action='store_true', help='Write residual file.')
    parser.add_argument(
        '--in_cluster', action='store_true', help='launch as a job in cluster.')

    ops = parser.parse_args()
    if not ops.relax_ic_in_csf:
        if ops.relax_ic_coeff is not None:
            raise ValueError(
                'Cannot use relax_ic_coeff when not using --relax_ic_in_csf')
    if ops.use_dti:
        if ops.in_dti is None:
            ops.in_dti = raw_input('Enter a DTI file (since you are using '
                                   '--use_dti.):  ')
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
    with open(ops.petsc_op_file, 'r') as fil:
        # Get all lines in the file that does not start with '#'
        # into a list
        petsc_ops = [
            line.strip() for line in fil if not line.startswith('#')]

    optional_args, bool_args = [], []
    if ops.div12pt_stencil:
        bool_args.append('--div12pt_stencil')
    if ops.relax_ic_in_csf:
        bool_args.append('--relax_ic_in_csf')
        if ops.relax_ic_coeff is not None:
            optional_args.append('-relax_ic_coeff ' + ops.relax_ic_coeff)
    if ops.mu_file is not None:
        optional_args.append('-muFile ' + ops.mu_file)
    if ops.use_dti:
        bool_args.append('--useTensorLambda')
        optional_args.append('-lambdaFile ' + op.join(res_dir, ops.in_dti))
    if ops.wrt_press:
        bool_args.append('--writePressure')
    if ops.wrt_force:
        bool_args.append('--writeForce')
    if ops.wrt_residual:
        bool_args.append('--writeResidual')

    cmd = ('%s -parameters %s -boundary_condition %s -atrophyFile %s '
           '-maskFile %s -imageFile %s -numOfTimeSteps %s -resPath %s '
           '-resultsFilenamesPrefix %s %s'
           % (target, ops.lame_paras, ops.boundary_condition, atrophy,
              in_seg, in_img, ops.time_steps, res_path, ops.res_prefix,
              ' '.join(bool_args + optional_args + petsc_ops)))

    if ops.in_cluster is True:
        cluster_mpi = '/opt/openmpi-gcc/current/bin/mpiexec '
        cmd = bu.sophia_nef_pbs_setting() + cluster_mpi + cmd
        job_name = '%sSteps%s' % (ops.res_prefix, ops.time_steps)
        #print cmd
        bu.qsub_job(name=job_name, queue='parlong', walltime='12:00:00',
                    procs='nodes=3:xeon:ppn=20', mem='mem=540gb',
                    #procs='nodes=4:xeon:ppn=20', mem='mem=720gb',
                    #procs='nodes=5:xeon:ppn=20', mem='mem=900gb',
                    dest_dir=res_dir, cmd=cmd)
    else:
        #print cmd
        bu.print_and_execute(cmd)


if __name__ == "__main__":
    main()
