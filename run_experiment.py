#!/usr/bin/env python
"""
file that ultimately guides a model through the simulation 
and starts the postprocessing
"""
from accessories import *
from fenics import *
from datetime import datetime
from functools import partial
import numpy as np
import threading
import os


def run_experiment(exp, mod, submodel, path="../", path_bias="", silent=False, **kwargs):
    # filter out kwargs which are needed in this method
    if "save_freq" in kwargs: 
        save_freq=kwargs.pop("save_freq")
    else: save_freq = -1.0
    if "timectrl" in kwargs:
        timectrl = kwargs.pop("timectrl")
    else: timectrl= False
    if "bot" in kwargs:
        bot = kwargs.pop("bot")
    else:
        bot=None

    # setting save_freq to -1.0 will later always fulfill the if condition
    model = mod(submodel, silent = silent)
    experiment = exp(**kwargs)
    model.init_experiment(experiment)
    model.silent = silent
    
    
    timestamp = datetime.now()
    timestamp = str(timestamp.day)+"-"+str(timestamp.month)+"-"+str(timestamp.year)+"_"+str(timestamp.hour)+"-"+str(timestamp.minute)
    if not silent: print("The current timestamp is: ",timestamp)
    if bot != None: bot.send_message("Started simulation with timestamp "+timestamp)
    path = path+"simulations/"+path_bias+"_sim_"+timestamp
    os.makedirs(path,exist_ok=True)


    parameters["form_compiler"]["precision"] = 100

    # - init logs
    logs = log_class(path)
    param_log = logs.add_log("parameter_log")
    time_step_log =logs.add_log("time_scheme_log", sort_by=["time"])
    dim = experiment.dim
    time_arr = []
    Energies=[]
    # - how many decimal points to be rounded for saving
    decimal_places = int(np.ceil(-1*np.log10(model.dt)))


    # - logging parameters
    param_log.add({**{"model name":model.modelname, "experiment": experiment.name ,"timestamp":timestamp, "maximum iterations" : "TODO", "tolerance of fixpoint solver": "TODO", "T":model.T,"dt":model.dt,"dh":model.dh, "plot frequency":save_freq, "time step control":timectrl},**experiment.param_dict})
    param_log.save(path)
    # - init postprocess 
    postprocess_v = fem_postprocess( dim, ["xdmf","quiver"], path, "v")
    postprocess_v.quiver_scale = (0.0,2.0)
    postprocess_v.decimal_places = decimal_places
    postprocess_d = fem_postprocess( dim, ["xdmf","quiver"], path, "d")
    postprocess_d.quiver_scale = (-1.0,1.0)
    postprocess_d.abs = False
    postprocess_d.decimal_places = decimal_places
    # - save initial conditions
    (vl,pl,dl,ql)=model.get_functions()
    postprocess_v.save(vl,experiment.mesh, 0.0)
    postprocess_d.save(dl,experiment.mesh, 0.0)

    if not silent: print("Starting the time evolution")
    model.computation_time = 0.0
    try:
        while model.t <= model.T:
            model.t += model.dt
            if not silent: print("Starting time-step ", model.t)
            (iteration, fp_errs)=model.do_time_step()
            # evtl. use time step control
            if timectrl: model.update_time_step_size()
            # set solution of current iteration as initial condition of next iteration
            model.update_ics()
            
            # - logging errors
            log_dict = {"time": model.t , "iterations": iteration ,  "fp_err":dict(zip(["v","d","q"], fp_errs)) , "processing_time": model.computation_time }
            # the following can be activated in order to save certain metrics to the log
            log_dict = { **log_dict, **model.get_model_dependent_log()}
            time_step_log.add(log_dict)
            time_step_log.save(path)

            # - postprocessing
            if model.t< ((model.dt)*3/2) or ((model.t+model.dt)%save_freq) < (3*model.dt/2):
                # the following sends bot.photo_freq photos via the telegram channel
                if bot!=None and (np.abs(model.t % (model.T/bot.photo_freq))< np.maximum(save_freq,model.dt)):
                    send_bot = bot
                else: 
                    send_bot = None
                postprocess_v.save(vl,experiment.mesh, model.t, send_bot=send_bot)
                postprocess_d.save(dl,experiment.mesh, model.t, send_bot=send_bot)
                
                time_arr.append(model.t)
                Energies.append(model.energy)
                energy_plot(time_arr, Energies,model.energy_labels, path)
                if bot!=None and (np.abs(model.t % (model.T/bot.photo_freq))< np.maximum(save_freq,model.dt)): 
                    bot.send_photo(path+"/plots/energy_plot.png")
                    bot.send_message("Finished "+str(round(model.t/model.T*100,0))+" percent of the simulation with the timestamp "+timestamp)

                

            model.update_ics()
            
    except Exception as err:
        print("Warning: Ran into error", err)
        param_log.add({"Failed": err })
        param_log.save(path)
      
    logs.save_all()


def run_multiple_experiments(cluster, sim_list, path="../", no_of_threads=1, silent=False):
    """
    The run_id gives all simulations the same characters in the beginning of the folder name for easier comparison.
    """
    if not silent: print("no of threads selected for parallel simulations: ", no_of_threads)
    tasks=[]
    threadlist=[]
    editing_thread = 0
    run_id = 1
    for (exp, mod, submodel, kwargs) in sim_list:
        
        # if not silent: print("loading experiment...")
        # exp = e(dim=dim, dt=dt, n=dh)
        # if test_run == True: 
        #     exp.T=2*dt
        # if not silent: print("initializing experiment setting...", exp.name, " for model ", model.modelname)
        # model.init_experiment(exp)
        # s="running: "+", ".join([model.modelname, exp.name,"dt", str(dt), "dh", str(dh)])
        # print(s)
        t=partial(run_experiment, exp, mod, submodel,path=path, path_bias=str(cluster+"_run"+str(run_id)), silent = silent, **kwargs)
        
        threadlist.append(threading.Thread(target=t))
        print("task to threadlist added...")
        # mod no_of_threads
        print("current thread ", editing_thread)
        if editing_thread==(no_of_threads-1):
            for threadinstance in threadlist:
                print("starting thread...")
                threadinstance.start()
            for threadinstance in threadlist:
                print("waiting for thread...")
                threadinstance.join()
            threadlist=[]
            editing_thread=-1
            try: 
                summary(cluster)
            except:
                print("Warning: summary can currently not be accessed...")
        tasks.append((t,"test"))
        run_id+=1
        editing_thread+=1