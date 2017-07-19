#!/usr/bin/env python 
#  -*- coding:utf-8 -*-
################################################################
# MODULE: r.massmov.sensitivity 
# AUTHOR(S): Monia Molinari, Massimiliano Cannata
# PURPOSE: sensitivity analysis tool for r.massmov model
# COPYRIGHT: (C) 2012 by the GRASS Development Team 
# 
# This program is free software under the GNU General Public 
# License (>=v2). Read the file COPYING that comes with GRASS 
# for details. 
# ##############################################################
#%Module
#%  description: sensitivity analysis tool for r.massmov model
#%  keywords: raster, landslide, modeling
#%End

#%option 
#% key: elev 
#% type: string 
#% gisprompt: old,cell,raster 
#% description: Name of elevation raster map
#% required: yes 
#%end

#%option 
#% key: hini
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of landslide initial body thickness raster map
#% required: yes 
#%end

#%option 
#% key: fluiddist
#% type: string 
#% gisprompt: old,cell,raster 
#% description: Distance to the toe of the landslide raster map
#% required: yes
#%end

#%option 
#% key: ifrict 
#% type: double 
#% description: Angle of internal friction [degrees]
#% required: yes
#% multiple: yes 
#%end

#%option 
#% key: fluid 
#% type: double 
#% description: Upward velocity of transition from solid to fluid of the landsliding mass [m/s]
#% required: yes
#% multiple: yes 
#%end

#%option 
#% key: chezy 
#% type: double 
#% description: Chezy roughness coefficient [m/s2]
#% required: no
#% multiple: yes
#% guisection: In_optional
#%end

#%option 
#% key: bfrict 
#% type: double 
#% description: Angle of basal friction [degrees]
#% required: no
#% multiple: yes
#% guisection: In_optional
#%end

#%option 
#% key: timesteps
#% type: integer 
#% description: Maximum number of time steps of the simulation [s]
#% required: yes
#%end

#%option 
#% key: stop_thres
#% type: double 
#% description: Pearson value threshold for simulation stop [-]
#% required: no
#% guisection: In_optional
#%end

#%option 
#% key: step_thres
#% type: integer 
#% description:  Number of time steps for evaluating stop_thres value [-]
#% required: no
#% guisection: In_optional
#%end

#%option 
#% key: threads
#% type: integer 
#% description: Number of threads for parallel computing [-]
#% required: no
#% guisection: In_optional
#%end

#%option
#% key: color_rules
#% type: string
#% gisprompt: new_file,file,input
#% description: Path to color rules file
#% required: no
#% guisection: In_optional
#%end

#%option 
#% key: prefix
#% type: string 
#% gisprompt: new,cell,raster 
#% description: Prefix for flow thickness output raster maps
#% required: no
#% guisection: Out_Optional
#%end

#%option
#% key: outfile
#% type: string
#% gisprompt: new_file,file,input
#% description: Path to report output file
#% required: yes
#%end

#%flag
#% key: m
#% description: Calculate mean error between simulation and observation
#% guisection: Out_Optional
#%end

#%flag
#% key: s
#% description: Calculate standard deviation between simulation and observation
#% guisection: Out_Optional
#%end

#%flag
#% key: c
#% description: Calculate Pearson correlation coefficient
#% guisection: Out_Optional
#%end

#%flag
#% key: p
#% description: Plot profile sets in xy graph 
#% guisection: Out_Optional
#%end

#%option 
#% key: profile
#% type: string 
#% description: Semicolon separated profiles given by one or more coordinate pairs (P11x,P11y[[,P12x,P12y,...];P21x,P21y[,P22x,P22y,...]])
#% required: no
#% multiple: yes
#% guisection: In_optional
#%end

#%option 
#% key: out_profile
#% type: string 
#% gisprompt: new_file,file,input 
#% description: Path to profile reports 
#% required: no
#% guisection: Out_Optional
#%end

import re
import sys
import time
import math
import numpy
import tempfile
import csv
import pylab
import matplotlib.pyplot as plt
import grass.script as grass


def ExtList(li,med):
    if med in li:
        pass
    else:
        li.append(med)
    return sorted(li)


def GetInput(par,name_par):    
    if (par.find(":")==-1 and par.find(",")==-1):
        try:
            isinstance(float(par),float)
            list_r=[]
            list_r.append(float(par))
            return list_r
        except:
            pass
           
    elif (par.count(":")==2):            
        list_par = par.split(":")
        for item in list_par:
            try:
                isinstance(float(item),float)
            except:
                 grass.fatal(_("%s parameter input format uncorrect")%name_par )
        list_r=[]
        r = float(list_par[0])
        while r <= float(list_par[1]):
            list_r.append(r)
            r+=float(list_par[2])
        if list_r==[]:
            grass.fatal(_("No values specified for %s parameter")%name_par )
        return list_r             
    elif (par.find(",")!=-1):
        list_par = par.split(",")
        list_ok=[]
        for item in list_par:
            try:
                isinstance(float(item),float)
                list_ok.append(float(item))
            except:
                 grass.fatal(_("%s parameter input format uncorrect")%name_par )   
        return list_ok   
    else:
        grass.fatal(_("%s parameter input format uncorrect")%name_par )    


def GetProfileCoord(profile):
    list_profile=[]
    if (profile.find(";")>=1):
        x = str(profile).split(";")
        for item in range(len(x)):
            list_profile.append(x[item])
        return list_profile
    elif (profile.find(",")>=4):
        list_profile.append(profile)
        return list_profile
    else:
        grass.fatal(_("Profiles coordinate values input format uncorrect") )    
    

def GetPlot(p_set,p_num,list_par,name_par,out_path):
    num_plots=len(list_par)
    for x in range(p_num):
        k=0        
        fig = plt.figure(figsize=(12,6), dpi=300)        
        # Set colormap
        colormap = plt.cm.jet
        plt.gca().set_color_cycle([colormap(i) for i in numpy.linspace(0.0, 0.9, num_plots)])                
        ax = plt.subplot(111)              
        for i in range (x,len(p_set),p_num):
            ax.plot(p_set[i][0],p_set[i][1],label='%s'%list_par[k])
            k+=1        
        # Set label and axes   
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel('Deposit height (m)')
        ax.set_title('Profile %s - %s'%(x,name_par))
        ax.grid(True)
        # Shink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        # Set legend       
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),fancybox=True, shadow=True)    
        pylab.savefig('%sProfile%s-%s'%(out_path,x,name_par))
        pylab.clf()
        pylab.cla()
     
def GetPearson(obs,model):
    p=re.compile('\W')
    num = p.sub("",tempfile.mktemp())
    den1 = p.sub("",tempfile.mktemp())    
    den2 = p.sub("",tempfile.mktemp()) 
    # Mean obs map #
    stats = grass.read_command("r.univar", map="%s" %obs,quiet=True)     
    mean_line_obs = stats.split("\n")[9]
    mean_obs = mean_line_obs.split(":")[1]    
    # Mean model map #
    stats = grass.read_command("r.univar", map="%s" %model,quiet=True)     
    mean_line_mod = stats.split("\n")[9]
    mean_mod = mean_line_mod.split(":")[1]    
    # Num (x1-xm)(y1-ym) #
    grass.mapcalc("$out1 = (($rast1-$m1)*($rast2-$m2))",out1=num,rast1=obs,rast2=model,\
    m1=mean_obs,m2=mean_mod,quiet=True)
    stats = grass.read_command("r.univar", map="%s" %num,quiet=True)
    sum_line_num = stats.split("\n")[14]
    summa_num = sum_line_num.split(":")[1]
    # Den1 sqrt(x1-xm)2 #   
    grass.mapcalc("$out2 = ($rast1-$m1)^2",out2=den1,rast1=obs,m1=mean_obs,quiet=True)
    stats = grass.read_command("r.univar", map="%s" %den1,quiet=True)
    sum_line_den = stats.split("\n")[14]
    summa_den1 = sum_line_den.split(":")[1]   
    # Den2 sqrt(y1-ym)2#
    grass.mapcalc("$out3 = ($rast2-$m2)^2",out3=den2,rast2=model,m2=mean_mod,quiet=True)
    stats = grass.read_command("r.univar", map="%s" %den2,quiet=True)
    sum_line_den = stats.split("\n")[14]
    summa_den2 = sum_line_den.split(":")[1]    
    p = float(summa_num)/(math.sqrt(float(summa_den1))*math.sqrt(float(summa_den2)))
    grass.run_command("g.remove",rast="%s,%s,%s"%(num,den1,den2),quiet=True)    
    return p

def GetMeanError(obs,model): 
    p=re.compile('\W')   
    diff = p.sub("",tempfile.mktemp())
    grass.mapcalc("$x = if(($obs!=0 || $mod!=0),($obs-$mod),null())",overwrite=True,x=diff,\
    obs=obs,mod=model,quiet=True)
    stats = grass.read_command("r.univar", map="%s"%diff,quiet=True)
    mean_line = stats.split("\n")[10]
    mean = mean_line.split(":")[1]    
    grass.run_command("g.remove",rast="%s"%diff,quiet=True)    
    return mean

def GetSteps(map_info):
    info = grass.read_command("r.info", map="%s"%map_info,quiet=True)
    step_line = info.split("\n")[20]
    step_white = re.sub( '\s+', ' ', step_line ).strip()
    step = step_white.split(" ")[3]
    return step

def GetSQMError(obs,model):  
    p=re.compile('\W') 
    diff = p.sub("",tempfile.mktemp())
    grass.mapcalc("$x = if(($obs!=0 || $mod!=0),($obs-$mod),null())",overwrite=True,\
    x=diff,obs=obs,mod=model,quiet=True)
    stats = grass.read_command("r.univar", map="%s"%diff,quiet=True)
    sqm_line = stats.split("\n")[11]
    sqm = sqm_line.split(":")[1]
    grass.run_command("g.remove",rast="%s"%diff,quiet=True)    
    return sqm

def main():
    elev = options["elev"]
    hini = options["hini"]
    fluiddist = options["fluiddist"]    
    ifrict = options["ifrict"]
    fluid = options["fluid"]    
    chezy = options["chezy"]   
    bfrict = options["bfrict"]  
    timesteps = int(options["timesteps"])  
    stop_thres = options["stop_thres"]
    step_thres = options["step_thres"]   
    threads =  options["threads"]
    outfile = options["outfile"]
    h = options["prefix"]
    profile = options["profile"]
    out_profile = options["out_profile"]
    color_rules = options["color_rules"]
    rheology="Voellmy"

    #### Check Input data ####
    
    if(rheology=='Voellmy'):
        if(len(chezy)==0 or len(bfrict)==0):
            grass.fatal(_("For the selected rheology Chezy and \
            basal friction parameters are required") )
            
    #### Check plotting profiles input data ####   
    if flags ['p']:
        if (len(profile)==0 or len(out_profile)==0):
            grass.fatal(_("For the profiles plotting <profile> \
            and <out_profile> parameters are required") )
            
    #### Check profiles input data ####   
    if (((len(profile)==0) and (len(out_profile)>0)) \
    or((len(out_profile)==0)and (len(profile)>0))):
        grass.fatal(_("For the profiles output <profile> and \
        <out_profile> parameters are required") )        

    #### Get general parameters values ####
    list_fluid = GetInput(fluid,"fluid")
    med_fluid = numpy.mean([list_fluid[0],\
    list_fluid[(len(list_fluid)-1)]])
    li_fluid = ExtList(list_fluid,med_fluid)
    list_ifrict = GetInput(ifrict,"ifrict")
    med_ifrict = numpy.mean([list_ifrict[0],\
    list_ifrict[(len(list_ifrict)-1)]])
    li_ifrict = ExtList(list_ifrict,med_ifrict)
    
    #### VOELLMY: Get parameters values ####        
    if(rheology=='Voellmy'):
        list_chezy = GetInput(chezy,"chezy")
        med_chezy = numpy.mean([list_chezy[0],\
        list_chezy[(len(list_chezy)-1)]])
        li_chezy = ExtList(list_chezy,med_chezy)        
        list_bfrict = GetInput(bfrict,"bfrict")       
        med_bfrict = numpy.mean([list_bfrict[0],\
        list_bfrict[(len(list_bfrict)-1)]])
        li_bfrict = ExtList(list_bfrict,med_bfrict)

    ### Define temporary maps ###
    p=re.compile('\W') 
    h_name = p.sub("",tempfile.mktemp()) 
    obs = p.sub("",tempfile.mktemp())
    

    ### Create mean observation map ###
    if (len(stop_thres)>0):
        if (len(threads)==0):
            grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
            fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
            fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
            stop_thres="%s" %stop_thres, step_thres="%s" %step_thres, timesteps="%s"\
            %timesteps, h="%s" %obs,quiet=True,overwrite=True)
            grass.run_command("r.null",map="%s"%obs, null=0)
        else:
            grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
            fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
            fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
            stop_thres="%s" %stop_thres, step_thres="%s" %step_thres, timesteps="%s"\
            %timesteps, h="%s" %obs,quiet=True,overwrite=True,threads="%s"%threads)
            grass.run_command("r.null",map="%s"%obs, null=0)            
    else:
        if (len(threads)==0):        
            grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
            fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
            fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
            timesteps="%s" %timesteps, h="%s" %obs,quiet=True,overwrite=True)
            grass.run_command("r.null",map="%s"%obs, null=0) 
        else:
            grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
            fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
            fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
            timesteps="%s" %timesteps, h="%s" %obs,quiet=True,overwrite=True,threads="%s"%threads)
            grass.run_command("r.null",map="%s"%obs, null=0)             
           

    ### Setting output file ###
    f = open("%s"%outfile,"w")
    f.write ("N_SIM"+","+"I_FRICT"+","+"FLUID_RATE"+","+"B_FRICT"+","+\
    "CHEZY"+","+"TIME"+","+"TIMESTEPS")
    if flags ['m']:
        f.write (","+"MEAN_ERR")
    if flags ['s']:
        f.write (","+"SQM_ERR")
    if flags ['c']:
        f.write (","+"PEARSON")
    f.write ("\n")

    count = 1

    #### VOELLMY: Launch combinations parameters values ####  
   
    if(rheology=='Voellmy'):   
        import time
        
        set_p=[]                
        #### IFRICT ####    
        for i in li_ifrict:
            if (len(stop_thres)>0):
                if (len(threads)==0):
                    start_time = time.time()
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %i,\
                    fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
                    stop_thres="%s" %stop_thres, step_thres="%s" %step_thres,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True)
                    elapsed = time.time() - start_time
                else:
                    start_time = time.time()
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %i,\
                    fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
                    stop_thres="%s" %stop_thres, step_thres="%s" %step_thres,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True, threads="%s"%threads)
                    elapsed = time.time() - start_time
            else:
                if (len(threads)==0):
                    start_time = time.time()
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %i,\
                    fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True) 
                    elapsed = time.time() - start_time 
                else: 
                    start_time = time.time()
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %i,\
                    fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True,threads="%s"%threads) 
                    elapsed = time.time() - start_time   
                               
            tot_steps = GetSteps(h_name)                        
            grass.run_command("r.null",map="%s"%h_name,null=0,quiet=True)
            if flags ['m']:
                MEAN = GetMeanError(obs,h_name)
            if flags ['s']:
                SQM = GetSQMError(obs,h_name)
            if flags ['c']:
                PEARS = GetPearson(obs,h_name)           
            f.write ("%s"%count+","+"%s"%i+","+"%s"%med_fluid+","+"%s"%med_bfrict+","\
            +"%s"%med_chezy+","+"%s"%elapsed+","+"%s"%tot_steps)
            if flags ['m']:
                f.write (","+"%s"%MEAN)
            if flags ['s']:
                f.write (","+"%s"%SQM)
            if flags ['c']:
                f.write (","+"%s"%PEARS)
            f.write("\n")
                
        
            ### Get Profile ###
            if len(h)>0:
                name_p_map=h
            else:
                name_p_map="h"
                
            if profile:
                list_p = GetProfileCoord(profile)
                                 
                for n_profile,item in enumerate(list_p):
                  
                    if (grass.run_command("r.profile",quiet=True,flags="g",\
                    input="%s" %(h_name),output="%sp%s_%s_%s.txt" \
                    %(out_profile,n_profile,name_p_map,count),profile="%s" %item)!=0):
                        grass.fatal(_('r.profile failed'))
                    
                    if flags['p']:
                        csv_reader = csv.reader(open('%sp%s_%s_%s.txt'%(out_profile,\
                        n_profile,name_p_map,count)),delimiter=' ')
                        p=[]
                        x=[]
                        y=[]
                        for line in csv_reader:
                            x.append(line[2])
                            y.append(line[3])
                            p=[x,y]
                        set_p.append(p)

            ### rename or delete temporary maps ### 
            if len(h)>0:                
                grass.run_command("g.rename",rast="%s,%s_%s"%(h_name,h,count),\
                overwrite=True,quiet=True)
                if (len(color_rules)>0):
                    grass.run_command("r.colors",map="%s_%s"%(h,count),\
                    rules="%s"%color_rules,quiet=True)
            else:
                grass.run_command("g.remove",rast="%s"%h_name,quiet=True)

            count+=1
       
        ### Get Plot ###                       
        if flags['p']:
            GetPlot(set_p,len(list_p),li_ifrict,'Internal Friction',out_profile)
                
        set_p=[]         
        #### FLUID ####     

        for ff in li_fluid:
            if (len(stop_thres)>0):  
                if (len(threads)==0):
                    start_time = time.time()    
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %ff, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
                    stop_thres="%s" %stop_thres, step_thres="%s" %step_thres,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True) 
                    elapsed = time.time() - start_time 
                else:
                    start_time = time.time()    
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %ff, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
                    stop_thres="%s" %stop_thres, step_thres="%s" %step_thres,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True,threads="%s"%threads) 
                    elapsed = time.time() - start_time                     
            else:
                if (len(threads)==0):
                    start_time = time.time()    
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %ff, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True) 
                    elapsed = time.time() - start_time 
                else:
                    start_time = time.time()    
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %ff, bfrict="%s" %med_bfrict, chezy="%s" %med_chezy,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True,threads="%s"%threads) 
                    elapsed = time.time() - start_time                    

 
            tot_steps = GetSteps(h_name)            
            grass.run_command("r.null",map="%s"%h_name,null=0,quiet=True)
            if flags ['m']:
                MEAN = GetMeanError(obs,h_name)
            if flags ['s']:
                SQM = GetSQMError(obs,h_name)
            if flags ['c']:
                PEARS = GetPearson(obs,h_name)
            f.write ("%s"%count+","+"%s"%med_ifrict+","+"%s"%ff+","+"%s"%med_bfrict+\
            ","+"%s"%med_chezy+","+"%s"%elapsed+","+"%s"%tot_steps)
            if flags ['m']:
                f.write (","+"%s"%MEAN)
            if flags ['s']:
                f.write (","+"%s"%SQM)
            if flags ['c']:
                f.write (","+"%s"%PEARS)
            f.write("\n")


            ### Get Profile ###

            if len(h)>0:
                name_p_map=h
            else:
                name_p_map="h"
                

            if profile:
                list_p = GetProfileCoord(profile)

                for n_profile,item in enumerate(list_p):
                    if (grass.run_command("r.profile",quiet=True,flags="g",\
                    input="%s" %(h_name),output="%sp%s_%s_%s.txt" \
                    %(out_profile,n_profile,name_p_map,count),profile="%s" %item)!=0):
                        grass.fatal(_('r.profile failed'))
                        
                    if flags['p']:
                        csv_reader = csv.reader(open('%sp%s_%s_%s.txt'%(out_profile,\
                        n_profile,name_p_map,count)),delimiter=' ')
                        p=[]
                        x=[]
                        y=[]
                        for line in csv_reader:
                            x.append(line[2])
                            y.append(line[3])
                            p=[x,y]
                        set_p.append(p)

            ### rename or delete temporary maps ### 
            if len(h)>0:
                grass.run_command("g.rename",rast="%s,%s_%s"%(h_name,h,count),\
                overwrite=True,quiet=True)
                if (len(color_rules)>0):
                    grass.run_command("r.colors",map="%s_%s"%(h,count),\
                    rules="%s"%color_rules,quiet=True)
            else:
                grass.run_command("g.remove",rast="%s"%h_name,quiet=True)
                        
            count+=1

        ### Get Plot ###  
        if flags['p']:
            GetPlot(set_p,len(list_p),li_fluid,'Fluid Rate',out_profile)
        
        set_p=[]         
        #### CHEZY ####  
        
        for c in li_chezy:
            if (len(stop_thres)>0):   
                if (len(threads)==0): 
                    start_time = time.time()           
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %c,\
                    stop_thres="%s" %stop_thres, step_thres="%s" %step_thres,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True) 
                    elapsed = time.time() - start_time 
                else:
                    start_time = time.time()           
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %c,\
                    stop_thres="%s" %stop_thres, step_thres="%s" %step_thres,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True,threads="%s"%threads) 
                    elapsed = time.time() - start_time                      
            else:
                if (len(threads)==0): 
                    start_time = time.time()           
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %c,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True) 
                    elapsed = time.time() - start_time 
                else:
                    start_time = time.time()           
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %med_fluid, bfrict="%s" %med_bfrict, chezy="%s" %c,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True,threads="%s"%threads) 
                    elapsed = time.time() - start_time                   


            tot_steps = GetSteps(h_name)           
            grass.run_command("r.null",map="%s"%h_name,null=0,quiet=True)
            if flags ['m']:
                MEAN = GetMeanError(obs,h_name)
            if flags ['s']:
                SQM = GetSQMError(obs,h_name)
            if flags ['c']:
                PEARS = GetPearson(obs,h_name)
            f.write ("%s"%count+","+"%s"%med_ifrict+","+"%s"%med_fluid+","+"%s"%med_bfrict+\
            ","+"%s"%c+","+"%s"%elapsed+","+"%s"%tot_steps)
            if flags ['m']:
                f.write (","+"%s"%MEAN)
            if flags ['s']:
                f.write (","+"%s"%SQM)
            if flags ['c']:
                f.write (","+"%s"%PEARS)
            f.write("\n")


            ### Get Profile ###           
            if len(h)>0:
                name_p_map=h
            else:
                name_p_map="h"
                
            if profile:
                list_p = GetProfileCoord(profile)

                for n_profile,item in enumerate(list_p):
                    if (grass.run_command("r.profile",quiet=True,flags="g",\
                    input="%s" %(h_name),output="%sp%s_%s_%s.txt" \
                    %(out_profile,n_profile,name_p_map,count),profile="%s" %item)!=0):
                        grass.fatal(_('r.profile failed'))

                    if flags['p']:
                        csv_reader = csv.reader(open('%sp%s_%s_%s.txt'\
                        %(out_profile,n_profile,name_p_map,count)),delimiter=' ')
                        p=[]
                        x=[]
                        y=[]
                        for line in csv_reader:
                            x.append(line[2])
                            y.append(line[3])
                            p=[x,y]
                        set_p.append(p)
                        
        ### rename or delete temporary maps ###  
            if len(h)>0:
                grass.run_command("g.rename",rast="%s,%s_%s"%(h_name,h,count),\
                overwrite=True,quiet=True)
                if (len(color_rules)>0):
                    grass.run_command("r.colors",map="%s_%s"%(h,count),\
                    rules="%s"%color_rules,quiet=True)
            else:
                grass.run_command("g.remove",rast="%s"%h_name,quiet=True)

            count+=1
            
        ### Get Plot ###  
        if flags['p']:
            GetPlot(set_p,len(list_p),li_chezy,'Chezy Coefficient',out_profile)

        set_p=[] 
        #### BFRICT ####  

        for b in li_bfrict:
            if (len(stop_thres)>0): 
                if(len(threads)==0):
                    start_time = time.time()    
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %med_fluid, bfrict="%s" %b, chezy="%s" %med_chezy,\
                    stop_thres="%s" %stop_thres, step_thres="%s" %step_thres,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True)  
                    elapsed = time.time() - start_time
                else:
                    start_time = time.time()    
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %med_fluid, bfrict="%s" %b, chezy="%s" %med_chezy,\
                    stop_thres="%s" %stop_thres, step_thres="%s" %step_thres,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True,threads="%s"%threads)  
                    elapsed = time.time() - start_time
            else:
                if(len(threads)==0):
                    start_time = time.time() 
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %med_fluid, bfrict="%s" %b, chezy="%s" %med_chezy,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True)   
                    elapsed = time.time() - start_time 
                else:
                    start_time = time.time() 
                    grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini,\
                    fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %med_ifrict,\
                    fluid="%s" %med_fluid, bfrict="%s" %b, chezy="%s" %med_chezy,\
                    timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True,threads="%s"%threads)   
                    elapsed = time.time() - start_time                


            tot_steps = GetSteps(h_name)            
            grass.run_command("r.null",map="%s"%h_name,null=0,quiet=True)
            if flags ['m']:
                MEAN = GetMeanError(obs,h_name)
            if flags ['s']:
                SQM = GetSQMError(obs,h_name)
            if flags ['c']:
                PEARS = GetPearson(obs,h_name)
            f.write ("%s"%count+","+"%s"%med_ifrict+","+"%s"%med_fluid+","+"%s"%b+","\
            +"%s"%med_chezy+","+"%s"%elapsed+","+"%s"%tot_steps)
            if flags ['m']:
                f.write (","+"%s"%MEAN)
            if flags ['s']:
                f.write (","+"%s"%SQM)
            if flags ['c']:
                f.write (","+"%s"%PEARS)
            f.write("\n")


            ### Get Profile ###

            if len(h)>0:
                name_p_map=h
            else:
                name_p_map="h"
            



            if profile:
                list_p = GetProfileCoord(profile)
                for n_profile,item in enumerate(list_p):
                    if (grass.run_command("r.profile",quiet=True,flags="g",\
                    input="%s" %(h_name),output="%sp%s_%s_%s.txt" %(out_profile,n_profile,\
                    name_p_map,count),profile="%s" %item)!=0):
                        grass.fatal(_('r.profile failed'))

                    if flags['p']:
                        csv_reader = csv.reader(open('%sp%s_%s_%s.txt'%(out_profile,n_profile,\
                        name_p_map,count)),delimiter=' ')
                        p=[]
                        x=[]
                        y=[]
                        for line in csv_reader:
                            x.append(line[2])
                            y.append(line[3])
                            p=[x,y]
                        set_p.append(p)

            ### rename or delete temporary maps ###       

            if len(h)>0:
                grass.run_command("g.rename",rast="%s,%s_%s"%(h_name,h,count),\
                overwrite=True,quiet=True)
                if (len(color_rules)>0):
                    grass.run_command("r.colors",map="%s_%s"%(h,count),\
                    rules="%s"%color_rules,quiet=True)
            else:
                grass.run_command("g.remove",rast="%s"%h_name,quiet=True)

            count+=1
            
        ### Get Plot ###            
        if flags['p']:
            GetPlot(set_p,len(list_p),li_bfrict,'Basal Friction',out_profile)

    ### rename or delete reference deposit map ###       
    grass.run_command("g.remove",rast="%s"%obs,quiet=True)
    
if __name__ == "__main__":
    options,flags = grass.parser()
    sys.exit(main())
