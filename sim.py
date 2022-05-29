#!/usr/bin/env python
"""
Main file: main method runs a single simulation specified by command_line inputs. For details type the following command into the terminal:
python main.py -h
example usage:
python sim.py -m decoupled_fp_solver -e annihilation_2 -s simple -d 3 -dt 0.01
--------
TO-DO: 
option for parallel computing
kwargs for run experiment via command line
same for model kwargs
"""
from run_experiment import *
from accessories import *
from lib import *
from fenics import *
from datetime import datetime
import numpy as np
import warnings
from argparse import ArgumentParser

# turn off deprecation warning
from ffc.quadrature.deprecation import QuadratureRepresentationDeprecationWarning
warnings.simplefilter("ignore", QuadratureRepresentationDeprecationWarning)


def main():
    """
    Definition of command line input and call of the simulation
    """
    parser = ArgumentParser(description="This program runs a single simulation with specification by command line input.")

    # necessary arguments
    parser.add_argument('-m','--mod', type=str, required=True, metavar='', help='THIS ARGUMENT IS REQUIRED: The model considered in the simulation. Options are '+", ".join(get_model_dict().keys()))
    parser.add_argument('-e','--exp', type=str, required=True, metavar='', help='THIS ARGUMENT IS REQUIRED: The experiment considered in the simulation. Options are '+", ".join(get_exp_dict().keys()))
    # default arguments
    parser.add_argument('-s','--submod', type=str, metavar='', nargs='?', const="simple", default="simple", help='The ericksen-leslie model considered in the simulation. Choices are "simple" or "general". Default is "simple"')
    parser.add_argument('-f','--freq', type=float, metavar='', nargs='?', const=-1.0, default=-1.0, help='Specifies the frequency with which the postprocessing is applied. Default is after every time step, i.e. freq=dt.')
    # optional arguments
    parser.add_argument('-d','--dim', type=int, metavar='', help='Configures the dimension of the experiment, if dim is supported by the chosen experiment')
    parser.add_argument('-dt','--dt', type=float, metavar='', help='Specifies resolution of time partition. Default is given by the experiment.')
    parser.add_argument('-dh','--dh', type=int, metavar='', help='Specifies resolution of space partition. Default is given by the experiment.')
    parser.add_argument('-i', '--id', type=str, metavar='', help='String that uniquely identifies this simulation.')
    parser.add_argument('-tc', '--timectrl', type=bool, metavar='', help='Specifies if a time step control shall be applied for the model.')
    #parser.add_argument('-k', '--kwargs', type=str, metavar='', help='Further optional arguments for the model.')
    parser.add_argument('-q', '--quiet', type=bool, metavar='', help='Turns off information in the terminal. Is of type boolean.')
    parser.add_argument('-p', '--path', type=str, metavar='', help='The path to the location where the simulation results will be saved. Default is the folder of main.py.')
    parser.add_argument('-not', "--notification",type=bool, metavar='', help='Sends out Telegram Notifications based on the file telebot_config.json.')
    # parse arguments
    args = parser.parse_args()
    # collect kwargs
    kwargs = dict()
    if args.dim: kwargs["dim"] = args.dim
    if args.dt: kwargs["dt"] = args.dt
    if args.dh: kwargs["n"] = args.dh
    if args.freq: kwargs["save_freq"] = args.freq
    if args.path: kwargs["path"] = args.path
    if args.id: kwargs["path_bias"] = args.id
    if args.quiet: kwargs["silent"] = args.quiet
    if args.timectrl: kwargs["timectrl"] = args.timectrl
    if args.notification and args.notification == True: 
        bot = msg_bot("telebot_config.json")
        kwargs["bot"] = bot

    model_dict = get_model_dict()
    exp_dict = get_exp_dict()
    # start simulation
    run_experiment(exp_dict[args.exp], model_dict[args.mod], args.submod, **kwargs)
    # optionally send stuff via telegram
    if args.notification and args.notification == True: bot.send_message("Simulation finished.")
        

    
    
if __name__ == "__main__":
    main()