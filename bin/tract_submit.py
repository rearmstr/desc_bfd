#! /usr/bin/env python

import time,os,argparse
from collections import deque
import subprocess as sub
import re
import os

def idSplit(id):
    if id is None:
        return id
    ids = []
    for r in id.split("^"):
        m = re.match(r"^(\d+)\.\.(\d+):?(\d+)?$", r)
        if m:
            limits = [int(v) if v else 1 for v in m.groups()]
            limits[1] += 1
            ids += range(*limits)
        else:
            ids.append(int(r))
    return ids
master=[8278,8279,8280,8281,8282,8283,8284,8285,8286,8519,8520,8521,8522,8523,8524,
        8525,8526,8527,8761,8762,8763,8764,8765,8766,8767,8768,8769,9003,9004,9005,
        9006,9007,9008,9009,9010,9070,9071,9072,9073,9074,9075,9076,9077,9078,9079,
        9102,9103,9104,9105,9106,9107,9125,9126,9127,9128,9129,9130,9131,9132,9133,
        9134,9135,9206,9207,9208,9209,9210,9211,9212,9213,9214,9215,9246,9247,9248,
        9249,9312,9313,9314,9315,9316,9317,9318,9319,9320,9321,9322,9344,9345,9346,
        9347,9348,9349,9350,9368,9369,9370,9371,9372,9373,9374,9375,9376,9377,9378,
        9448,9449,9450,9451,9452,9453,9454,9455,9456,9457,9458,9459,9555,9556,9557,
        9558,9559,9560,9561,9562,9563,9564,9565,9566,9587,9588,9589,9590,9591,9592,
        9593,9611,9612,9613,9614,9615,9616,9617,9618,9619,9620,9621,9691,9692,9693,
        9694,9695,9696,9697,9698,9699,9700,9701,9702,9798,9799,9800,9801,9802,9803,
        9804,9805,9806,9807,9808,9830,9831,9832,9833,9834,9835,9855,9856,9857,9858,
        9859,9860,9861,9933,9934,9935,9936,9937,9938,9939,9940,9941,9942,9943,9944,
        10040,10041,10042,10043,10044,10045,10046,10047,10048,10049,10176,10177,10178,
        10179,10180,10181,10182,10183,10185,10287,10288,10289,10290,10291,15827,15828,
        15829,15830,15831,15832,15833,16005,16006,16007,16008,16009,16010,16011,16012,
        16182,16183,16184,16185]

xmm=[8038,8039,8278,8279,8280,8281,8282,8283,8284,8285,8286,8519,8520,8521,8522,8523,8524,8525,8526,8527,8761,8761,8762,8763,8764,8765,8766,8767,8768,8769,9003,9004,9005,9006,9007,9008,9009,9010,9011,9245,9246,9247,9248,9249,9250,9251,9252,9253,9488,9489,9490,9491,9492,9493,9494,9495,9732,9733,9734,9735,9736,9737]
gama9=[9069,9070,9071,9072,9073,9074,9075,9076,9077,9078,9079,9080,9081,9082,9083,9312,9313,9314,9315,9316,9317,9318,9319,9320,9321,9322,9323,9324,9325,9326,9555,9556,9557,9558,9559,9560,9561,9562,9563,9564,9565,9566,9567,9568,9569,9797,9798,9799,9800,9801,9802,9803,9804,9805,9806,9807,9808,9809,9810,9811,10039,10040,10041,10042,10043,10044,10045,10046,10047,10048,10049,10282,10283,10284,10285,10286,10287,10288,10289,10290,10291,10292]
wide12=[9096,9097,9098,9099,9100,9101,9102,9103,9104,9105,9106,9107,9108,9109,9110,9111,9112,9113,9114,9338,9339,9340,9341,9342,9343,9344,9345,9346,9347,9348,9349,9350,9351,9352,9353,9354,9355,9356,9357,9581,9582,9583,9584,9585,9586,9587,9588,9589,9590,9591,9592,9593,9594,9595,9596,9597,9598,9599,9600,9581,9582,9583,9584,9585,9586,9587,9588,9589,9590,9591,9592,9593,9594,9824,9825,9826,9827,9828,9829,9830,9831,9832,9833,9834,9835,9836,9837,9838,9839,9840,9841,9842]
gama15=[9123,9124,9125,9126,9127,9128,9129,9130,9131,9132,9133,9134,9135,9136,9365,9366,9367,9368,9369,9370,9371,9372,9373,9374,9375,9376,9377,9378,9379,9611,9612,9613,9614,9615,9616,9617,9618,9619,9620,9621,9622,9851,9852,9853,9854,9855,9856,9857,9858,9859,9860,9861,9862,9863,9864]
hect=[15823,15824,15825,15826,15827,15828,15829,15830,15831,15832,15833,15834,16001,16002,16003,16004,16005,16006,16006,16007,16008,16009,16010,16011,16012,16176,16177,16178,16179,16180,16181,16182,16183,16184,16185,16185,16186]
vvds=[9206,9207,9208,9209,9210,9211,9212,9213,9214,9215,9216,9217,9218,9219,9448,9449,9450,9451,9452,9453,9454,9455,9456,9457,9458,9459,9460,9461,9462,9691,9692,9693,9694,9695,9696,9697,9698,9699,9700,9701,9702,9703,9704,9705,9933,9934,9935,9936,9937,9938,9939,9940,9941,9942,9943,9944,9945,9946,10176,10177,10178,10179,10180,10181,10182,10183,10184,10185,10186,10187,10188]

tract_dict = {
    'xmm': xmm,
    'gama9': gama9,
    'gama15': gama15,
    'hect': hect,
    'vvds': vvds,
    'wide12': wide12,
    'all': master
}

parser = argparse.ArgumentParser(description='Run single file')
parser.add_argument('--max_jobs',default=300,type=int,
                    help='number of jobs to run concurrently')
parser.add_argument('--njobs', default=50,type=int,help='number of jobs total to run')
parser.add_argument('--name',default='name',
                    help='name of submit job')
parser.add_argument('--exe',default=None,
                    help='executable')
parser.add_argument('--hours',type=int,default=8,
                    help='how many hours')
parser.add_argument('--mins',type=int,default=0,
                    help='how many minutes')
parser.add_argument('--mem',type=int,default=None,
                    help='memory needed in MB')
parser.add_argument('--n_threads',type=int,default=1,
                    help='number of threads')
parser.add_argument('--arg',default='',
                    help='generic arguments')
parser.add_argument('--tracts', default='8523',help='which tracts')
parser.add_argument('--field', default=None,
                               help='specify field from [all, xmm, hect, gama15, gama9, vvds, wide12]')
parser.add_argument('--bank', default=None,help='which bank to use')
parser.add_argument('--queue', default=None,help='which queue to use')

args = parser.parse_args()

if args.field is None:
    tracts = idSplit(args.tracts)
else:
    if args.field not in tract_dict.keys():
        raise('Not valid field: %s' % args.field)
    tracts = tract_dict[args.field]

tract_lists=[tracts[i::args.njobs] for i in range(args.njobs)]

# convert options to dictionary to use with format
dict=vars(args)
user = os.getlogin()

for ii,tract_list in enumerate(tract_lists):
    submit_list = '^'.join([str(tract) for tract in tract_list])
    while True:
        pipe = sub.Popen(['squeue','-u', '%s' % user ],stdout=sub.PIPE)
        # count the number of jobs currently running
        q_out=pipe.communicate()[0]
        num=len(str(q_out).split('\n'))-1
        if(num < args.max_jobs):break
        time.sleep(0.25)

    output='log.%d.%s'%(ii, args.name)
    dict['output']=output
    use_arg = dict['arg']
    use_arg += ' --id filter=HSC-I tract=%s' %(submit_list)

    dict['use_arg'] = use_arg

    dict['batch_bank'] = ''
    if args.bank is not None:
        dict['batch_bank'] = '#SBATCH -A %s' % args.bank

    dict['batch_queue'] = ''
    if args.queue is not None:
        dict['batch_queue'] = '#SBATCH -p %s' % args.queue

    dict['batch_mem'] = ''
    if args.mem is not None:
        dict['memory'] = '#SBATCH --mem %d' % args.mem

    submit_text="""#!/bin/bash
#SBATCH -N 1
#SBATCH -c {n_threads}
#SBATCH --output={output}
#SBATCH -t {hours}:{mins:02d}:00
{batch_mem}
{batch_bank}
{batch_queue}
{exe} {use_arg} """.format(**dict)


    submit_file='submit_%d_%s.cmd' % (ii,args.name)
    ofile=open(submit_file,'w')
    ofile.write(submit_text)
    ofile.close()
    os.system('sbatch %s' % (submit_file))



