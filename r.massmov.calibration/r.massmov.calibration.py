#!/usr/bin/env python
# -*- coding: utf-8 -*-
############################################################################ # 
# MODULE: r.massmov.calibration 
# AUTHOR(S): Monia Molinari, Massimiliano Cannata 
# PURPOSE: Uses r.mymodule to perform some elaboration 
# COPYRIGHT: (C) 2012 by the GRASS Development Team 
# 
# This program is free software under the GNU General Public 
# License (>=v2). Read the file COPYING that comes with GRASS 
# for details. 
# ############################################################################
#%Module
#%  description: Tool for the calibration of r.massmov model
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
#% key: obs 
#% type: string 
#% gisprompt: old,cell,raster 
#% description: Name of observation raster map
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
#% description: Number of time step for evaluating stop_thres value [-]
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

#%flag
#% key: c
#% description: Calculate Pearson correlation coefficient
#% guisection: Out_Optional
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


import re
import sys
import time
import numpy
import tempfile
import csv
import pylab
import math
import matplotlib.pyplot as plt
import grass.script as grass


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


def GetSteps(map_info):
    info = grass.read_command("r.info", map="%s"%map_info,quiet=True)
    step_line = info.split("\n")[20]
    step_white = re.sub( '\s+', ' ', step_line ).strip()
    step = step_white.split(" ")[3]
    return step

def GetMeanError(obs,model): 
    p=re.compile('\W')   
    diff = p.sub("",tempfile.mktemp())
    grass.mapcalc("$x = if(($obs!=0 || $mod!=0),($obs-$mod),null())",overwrite=True,x=diff,obs=obs,mod=model,quiet=True)
    stats = grass.read_command("r.univar", map="%s"%diff,quiet=True)
    mean_line = stats.split("\n")[10]
    mean = mean_line.split(":")[1]
    
    grass.run_command("g.remove",rast="%s"%diff,quiet=True)    
    return mean


def GetSQMError(obs,model):  
    p=re.compile('\W') 
    diff = p.sub("",tempfile.mktemp())
    grass.mapcalc("$x = if(($obs!=0 || $mod!=0),($obs-$mod),null())",overwrite=True,x=diff,obs=obs,mod=model,quiet=True)
    stats = grass.read_command("r.univar", map="%s"%diff,quiet=True)
    sqm_line = stats.split("\n")[11]
    sqm = sqm_line.split(":")[1]

    grass.run_command("g.remove",rast="%s"%diff,quiet=True)    
    return sqm

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
    grass.mapcalc("$out1 = (($rast1-$m1)*($rast2-$m2))",out1=num,rast1=obs,rast2=model,m1=mean_obs,m2=mean_mod,quiet=True)
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

def main():
    elev = options["elev"]
    obs = options["obs"]
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
    color_rules = options["color_rules"]
    profile = options["profile"]
    out_profile = options["out_profile"]
    rheology="Voellmy"



    #### Check Input data ####
    
    if(rheology=='Voellmy'):
        if(len(chezy)==0 or len(bfrict)==0):
            grass.fatal(_("For the selected rheology Chezy and basal friction parameters are required") )
      
   
    #### Get general parameters values ####
    
    li_fluid = GetInput(fluid,"fluid")
    
    li_ifrict = GetInput(ifrict,"ifrict")

    
    #### VOELLMY: Get parameters values #### 
    
    if(rheology=='Voellmy'):

        li_chezy=GetInput(chezy,"chezy")          
        li_bfrict=GetInput(bfrict,"bfrict") 



    ### Define temporary maps ###
    p=re.compile('\W') 
    h_name = p.sub("",tempfile.mktemp()) 


    ### Setting output file ###
    f = open("%s"%outfile,"w")
    f.write ("N_SIM"+","+"I_FRICT"+","+"FLUID_RATE"+","+"B_FRICT"+","+"CHEZY"+","+"TIME"+","+"TIMESTEPS")
    if flags ['m']:
        f.write (","+"MEAN_ERR")
    if flags ['s']:
        f.write (","+"SQM_ERR")
    if flags ['c']:
        f.write (","+"PEARSON")    
    f.write ("\n")

    count = 1

    # Get simulations foe each parameters combination #
    for xxx in li_ifrict:
        for yyy in li_bfrict:
            for zzz in li_chezy:
                for kkk in li_fluid:
                    if (len(stop_thres)>0):
                        if (len(threads)==0):
                            start_time = time.time()
                            grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini, fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %xxx,fluid="%s" %kkk, bfrict="%s" %yyy, chezy="%s" %zzz, stop_thres="%s" %stop_thres, step_thres="%s" %step_thres, timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True)
                            elapsed = time.time() - start_time  
                        else:
                            start_time = time.time()
                            grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini, fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %xxx,fluid="%s" %kkk, bfrict="%s" %yyy, chezy="%s" %zzz, stop_thres="%s" %stop_thres, step_thres="%s" %step_thres, timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True,threads="%s"%threads)
                            elapsed = time.time() - start_time 
                    else:
                        if (len(threads)==0):
                            start_time = time.time()
                            grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini, fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %xxx,fluid="%s" %kkk, bfrict="%s" %yyy, chezy="%s" %zzz, timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True)
                            elapsed = time.time() - start_time  
                        else:
                            start_time = time.time()     
                            grass.run_command ("r.massmov",elev="%s" %elev, h_ini="%s" %hini, fluiddist="%s" %fluiddist,rheology="Voellmy", ifrict="%s" %xxx,fluid="%s" %kkk, bfrict="%s" %yyy, chezy="%s" %zzz, timesteps="%s" %timesteps, h="%s" %h_name,quiet=True,overwrite=True,threads="%s"%threads)
                            elapsed = time.time() - start_time  


                    tot_steps = GetSteps(h_name) 
                    grass.run_command("r.null",map="%s"%h_name,null=0,quiet=True)
                    
                    if flags ['m']:
                        MEAN = GetMeanError(obs,h_name)
                    if flags ['s']:
                        SQM = GetSQMError(obs,h_name)
                    if flags ['c']:
                        P = GetPearson(obs,h_name)                        
                    

                    f.write ("%s"%count+","+"%s"%xxx+","+"%s"%kkk+","+"%s"%yyy+","+"%s"%zzz+","+"%s"%elapsed+","+"%s"%tot_steps)

                    if flags ['m']:
                        f.write (","+"%s"%MEAN)
                    if flags ['s']:
                        f.write (","+"%s"%SQM)
                    if flags ['c']:
                        f.write (","+"%s"%P)                    
                    f.write("\n")



                    ### Get profile with r.profile ###
                    if len(h)>0:
                        name_p_map=h
                    else:
                        name_p_map="h"
   
                    if profile:
                        list_p = GetProfileCoord(profile)
                                 
                        for n_profile,item in enumerate(list_p):
                  
                            if (grass.run_command("r.profile",quiet=True,flags="g",input="%s" %(h_name),output="%sp%s_%s_%s.txt" %(out_profile,n_profile,name_p_map,count),profile="%s" %item)!=0):
                                grass.fatal(_('r.profile failed'))

                    ### g.rename and r.colors on temporary map ### 
                    if len(h)>0:                
                        grass.run_command("g.rename",rast="%s,%s_%s"%(h_name,h,count),overwrite=True,quiet=True)
                        if (len(color_rules)>0):
                            grass.run_command("r.colors",map="%s_%s"%(h,count),rules="%s"%color_rules,quiet=True)

                    count+=1
    f.close()
      
    ## Delete temporary map ##
    if len(h)==0:
        grass.run_command("g.remove",rast="%s"%h_name,quiet=True)

if __name__ == "__main__":
    options,flags = grass.parser()
    sys.exit(main())
