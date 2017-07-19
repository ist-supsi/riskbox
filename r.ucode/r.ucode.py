#!/usr/bin/env python 
#  -*- coding:utf-8 -*-
############################################################################ # 
# MODULE: r.ucode
# AUTHOR(S): Monia Molinari, Massimiliano Cannata
# PURPOSE: UCODE-based sensitivity analysis and calibration for GRASS modules 
# COPYRIGHT: (C) 2012 by the GRASS Development Team 
# 
# This program is free software under the GNU General Public 
# License (>=v2). Read the file COPYING that comes with GRASS 
# for details. 
# ############################################################################
#%Module
#%  description: UCODE-based sensitivity analysis and calibration for GRASS modules
#%  keywords: raster, UCODE, calibration
#%End

#%option 
#% key: mode
#% type: string 
#% description: running mode
#% options: normal,config,advance
#% required: yes
#% answer: normal
#%end

#%option 
#% key: folder
#% type: string 
#% gisprompt: new_file,file,input 
#% description: Path to UCODE input/output files folder
#% required: yes
#%end

#%option 
#% key: obs_map 
#% type: string 
#% gisprompt: old,cell,raster 
#% description: Observation raster map
#% guisection: Obs_data
#% required: no 
#
#%end

#%option 
#% key: obs_points
#% type: string 
#% gisprompt: old,vector,vector 
#% description: Sampling points vector map
#% guisection: Obs_data
#% required: no 
#%end

#%option 
#% key: col_value
#% type: string 
#% description: column name of observations values
#% required: no
#% guisection: Obs_data
#%end

#%option 
#% key: col_name
#% type: string 
#% description: column name of observations names
#% required: no
#% guisection: Obs_data
#%end

#%option 
#% key: col_statflag
#% type: string 
#% description: column name of statistic used to calculate observations weight (VAR,SD,CV,WT,SQRWT)
#% required: no
#% guisection: Obs_data
#%end

#%option 
#% key: col_statistic
#% type: string 
#% description: column name of values used by selected statistic to calculate the observations weight
#% required: no
#% guisection: Obs_data
#%end

#%option 
#% key: module
#% type: string 
#% description: GRASS module name
#% required: no
#% guisection: Module_data
#%end

#%option 
#% key: fixed
#% type: string 
#% description: GRASS module fixed arguments
#% required: no
#% guisection: Module_data
#%end

#%option 
#% key: modout
#% type: string 
#% description: GRASS module raster output argument
#% required: no
#% guisection: Module_data
#%end

#%option 
#% key: par_table
#% type: string 
#% gisprompt: old_file,file,input 
#% description: GRASS module calibration parameters table 
#% guisection: Module_data
#% required: no
#%end

#%option 
#% key: prefix
#% type: string 
#% description: Prefix for ucode output files
#% required: no
#% guisection: Out_options
#%end

#%option 
#% key: exec
#% type: string 
#% gisprompt: old_file,file,input 
#% description: Path to ucode executable  
#% required: no
#% guisection: UCODE_options
#%end

#%option 
#% key: tolpar
#% type: double 
#% description: Tolerance based on parameter values
#% answer: 0.01
#% required: no
#% guisection: UCODE_options
#%end

#%option 
#% key: tolsosc
#% type: double 
#% description: Tolerance based on changes to model fit
#% answer: 0.0
#% required: no
#% guisection: UCODE_options
#%end

#%option 
#% key: maxiter
#% type: integer 
#% description: Maximum number of parameter-estimation iterations allowed
#% answer: 5
#% required: no
#% guisection: UCODE_options
#%end

#%option 
#% key: maxchange
#% type: double 
#% description: Maximum fractional amount parameter values allowed to change in one iteration
#% answer: 2.0
#% required: no
#% guisection: UCODE_options
#%end

#%option 
#% key: maxchange_realm
#% type: string 
#% description: Transformation space to which MaxChange applies
#% options: native,regression
#% required: no
#% answer: native
#% guisection: UCODE_options
#%end

#%flag
#% key: s
#% description: Sensitivity analysis mode
#% guisection: UCODE_options
#%end

#%flag
#% key: o
#% description: Calibration & sensitivity analysis mode
#% guisection: UCODE_options
#%end

#%flag
#% key: i
#% description: Omit Insensitive (only if -o activated)
#% guisection: UCODE_options
#%end


import re
import os
import sys
import time
import numpy
import tempfile
import fileinput
import csv
import pylab
import matplotlib.pyplot as plt
import grass.script as grass
import subprocess


def GetObsFile(rmap,points,col_n,col_v,col_s,col_sf):
    par = []
    attr = []
    max_len = []    
    

    ## ObsName column: extract name for each observation point ##
    if not len(col_n)>0:
        grass.run_command ("v.db.addcolumn", map="%s" %points, columns="ObsName varchar",quiet=True)
        grass.read_command("v.db.update", map="%s" %points, column="ObsName", qcolumn="'p_'||cat")
        list_name = grass.read_command("v.db.select", "c", map="%s" %points, columns="ObsName",quiet=True)
    else:
        list_name = grass.read_command("v.db.select", "c", map="%s" %points, columns="%s"%col_n,quiet=True)
       
    l_name = list_name.split("\n")[:-1]
    par.append(l_name)
    attr.append("obsname")
    max_len.append(len(max(l_name,key=len)))

    ## ObsValue column: extract observation points values ##    
    if not len(col_v)>0:
        grass.run_command ("v.db.addcolumn", map="%s" %points, columns="ObsValue DOUBLE",quiet=True)
        grass.run_command ("v.what.rast", map="%s" %points, column="ObsValue",raster="%s" %rmap,quiet=True)
        list_obsvalue = grass.read_command("v.db.select", "c", map="%s" %points, columns="ObsValue",quiet=True)
    else:
        list_obsvalue = grass.read_command("v.db.select", "c", map="%s" %points, columns="%s"%col_v,quiet=True)

    l_obsvalue = list_obsvalue.split("\n")[:-1]       
    par.append(l_obsvalue)
    attr.append("obsvalue")
    max_len.append(len(max(l_obsvalue,key=len)))

    ## Statist column: extract observation points weight ##   
    if len(col_s)>0:
        list_stat = grass.read_command("v.db.select", "c", map="%s" %points, columns="%s"%col_s,quiet=True)
        l_stat = list_stat.split("\n")[:-1]   
        par.append(l_stat)
        attr.append("statistic")
        max_len.append(len(max(l_stat,key=len)))
        
    ## StatFlag column: extract observation points statistic ##         
    if len(col_sf)>0:
        list_flag = grass.read_command("v.db.select", "c", map="%s" %points, columns="%s"%col_sf,quiet=True)
        l_flag = list_flag.split("\n")[:-1] 
        par.append(l_flag)
        attr.append("statflag")
        max_len.append(len(max(l_flag,key=len)))
            
    return (par,attr,max_len)

    
def GetInstrFile(points,proj_path,l_name):
    n_row = len(l_name)
    f = open("%sobs.instructions"%proj_path,"w")
    f.write("jif @\n")
    f.write("StandardFile 0 1 %s\n"%n_row)
    for item in l_name:
        f.write("%s\n"%item)
    f.close()
    
    
def GetTemplin(project,par_tpl):
    f = open("%stemplin.tpl"%project,"w")
    f.write("jtf @\n")
    for item in par_tpl:
        f.write("%s=@%s          @\n"%(item,item))
        
    f.close()
    
def main():
    obs_map = options["obs_map"]
    obs_points = options["obs_points"]
    module = options["module"]
    fixed = options["fixed"]
    modout = options["modout"]
    project = options["folder"]    
    par_table = options["par_table"]    
    exec_ucode = options["exec"] 
    prefix = options["prefix"] 
    tolsosc = float(options["tolsosc"]) 
    tolpar = float(options["tolpar"])  
    maxiter = int(options["maxiter"])  
    maxchange = float(options["maxchange"])     
    maxchange_realm = options["maxchange_realm"]  
    col_value =  options["col_value"]
    col_name =  options["col_name"] 
    col_statistic =  options["col_statistic"]   
    col_statflag =  options["col_statflag"]
    opti="yes"
    sens="no"
    omit="no"
    mode = options["mode"]
    

    ## set UCODE options ##              
    if flags ['s']:
        sens="yes"
        opti="no"
        
    if flags ['o']:
        sens="yes"
        
    if flags ['i']:
        omit="yes"


    ## check modalities ##
    #if (flags['g'] and flags['r']):
    #   grass.fatal(_('Only one mode <-g> or <-r> can be activated'))
        
    if (flags['s'] and flags['o']):
        grass.fatal(_('Only one mode <-o> or <-s> can be activated'))   
        
    if (flags['i'] and not flags['o']):
         grass.fatal(_('<-i> can be activated only for mode <-o>'))

    ## check required inputs ##
    
    #if not flags['r']:
    if not (mode=='advance'):
        if ((len(obs_points)>0) and (len(module)>0) and (len(fixed)>0) and (len(modout)>0) and (len(par_table)>0)):
            pass
        else:
            grass.fatal(_('For the selected mode parameters <obs_points>, <module>, <modout>, <fixed>,<par_table> are required'))
            
        if (((len(obs_points)>0)and (len(obs_map)>0))or((len(obs_points)>0)and (len(col_value)>0))):
            pass
        else:
            grass.fatal(_('For the selected mode parameters <col_value> or <obs_map> are needed'))
        
    #if flags['r']:
    if (mode=='advance'):
        if (len(obs_points)>0):
            pass
        else:
            grass.fatal(_('To obtain observations values one parameter between <col_value> and <obs_map> is required'))    
       
    ## normal o -g modes ##
    #if not flags['r']:
    if not (mode=='advance'):
        
        ## create sampling points vector map copy ##
        p=re.compile('\W') 
        points_copy = p.sub("",tempfile.mktemp())     
        grass.run_command ("g.copy", vect="%s,%s" %(obs_points,points_copy),quiet=True)  

        ## collect observation data ##
        (par,attr,max_len)=GetObsFile(obs_map,points_copy,col_name,col_value,col_statistic,col_statflag)
      
        ## get parameters to calibrate ##
        table = []
        for line in fileinput.input("%s"%par_table):
            table.append(line)
        row=len(table)-1
        col=len(table[0].split())
        first_col = []
        for item in table:
            first_col.append(item.split()[0])
        par_cal=first_col[1::]

        ## Create instruction.obs file ##    
        GetInstrFile(points_copy,project,par[0])

        ## Create templin.tpl file ##     
        GetTemplin(project,par_cal)
       
        ## Create module.sh file ##
        import os       
        m = os.open("%smodule.sh"%project,os.O_RDWR|os.O_CREAT)
        os.write(m,"eval `cat %smodin.sen`;\n"%project)
        os.write(m,"%s %s"%(module,fixed))
        for item in par_cal:
            os.write(m," %s=$%s"%(item,item))
        os.write(m," %s --o;\n"%modout)
        #if not flags['g']:
        if not (mode=='config'):
            os.write(m,"v.db.addcolumn map=%s column='UCODE_extv DOUBLE';\n"%points_copy)
            os.write(m,"r.null map=%s null=0;\n"%modout.split('=')[1])
            os.write(m,"v.what.rast map=%s column='UCODE_extv' raster=%s;\n"%(points_copy,modout.split('=')[1]))
            os.write(m,"v.db.select -c map=%s columns='UCODE_extv' > %smodout._os;\n"%(points_copy,project))
            os.write(m,"v.db.dropcolumn map=%s column='UCODE_extv';"%points_copy)
            
        else:
            os.write(m,"v.db.addcolumn map=%s column='UCODE_extv DOUBLE';\n"%obs_points)
            os.write(m,"r.null map=%s null=0;\n"%modout.split('=')[1])
            os.write(m,"v.what.rast map=%s column='UCODE_extv' raster=%s;\n"%(obs_points,modout.split('=')[1]))
            os.write(m,"v.db.select -c map=%s columns='UCODE_extv' > %smodout._os;\n"%(obs_points,project))
            os.write(m,"v.db.dropcolumn map=%s column='UCODE_extv';"%obs_points)        
        os.close(m)        

        ## create ucode.in file ##   
        i = open("%sucode.in"%project,"w")        
        separ='#'+'-'*40
        
        ## UCODE block ##     
        i.write("%s\n# UCODE \n%s\n\n"%(separ,separ))
        i.write("BEGIN Options	KEYWORDS\n")
        i.write("  Verbose=5\n")
        i.write("END Options\n\n")
        
        ## UCODE-CONTROL INFORMATION block ##   
        i.write("%s\n# UCODE-CONTROL INFORMATION \n%s\n\n"%(separ,separ))  
        i.write("BEGIN UCODE_CONTROL_DATA KEYWORDS\n")
        i.write("ModelName=%s\n"%module)
        i.write("sensitivities=%s\n"%sens)
        i.write("optimize=%s\n"%opti)
        i.write("OmitInsensitive=%s\n"%omit)
        i.write("END UCODE_CONTROL_DATA\n\n")
        i.write("BEGIN REG_GN_CONTROLS KEYWORDS\n")
        i.write("tolpar=%s\n"%tolpar)
        i.write("tolsosc=%s\n"%tolsosc)
        i.write("maxiter=%s\n"%maxiter)
        i.write("maxchange=%s\n"%maxchange)
        i.write("MaxchangeRealm=%s\n"%maxchange_realm)  
        i.write("END REG_GN_CONTROLS\n\n")
        
        ## COMMAND FOR APPLICATION MODEL(S) block ##   
        i.write("%s\n# COMMAND FOR APPLICATION MODEL(S) \n%s\n\n"%(separ,separ))  
        i.write("BEGIN MODEL_COMMAND_LINES TABLE\n")
        i.write("nrow=1 ncol=3  columnlabels\n")
        i.write("COMMAND    PURPOSE    COMMANDID\n")
        i.write("%smodule.sh    Forward       %s\n"%(project,module))
        i.write("END MODEL_COMMAND_LINES\n\n")       

        i.write("%s\n# PARAMETER INFORMATION \n%s\n\n"%(separ,separ))   
        i.write("BEGIN PARAMETER_DATA TABLE\n")
        i.write("nrow=%s ncol=%s columnlabels\n"%(row,col))
        for item in table:
            i.write("%s"%item)
        i.write("END PARAMETER_DATA\n\n")
        
        i.write("%s\n# OBSERVATION INFORMATION \n%s\n\n"%(separ,separ)) 
        i.write("BEGIN OBSERVATION_DATA TABLE\n")
        i.write("NROW=%s  NCOL=%s  COLUMNLABELS\n"%((len(par[0])),len(attr)))

        for item in attr:
            i.write("%s"%item.ljust(15))
        i.write("\n") 
        
        for x in range(len(par[0])):
            for y in range(len(attr)):
                i.write("%s"%par[y][x].ljust(max_len[y]+10))
            i.write("\n")
        i.write("END OBSERVATION_DATA\n\n") 

        i.write("%s\n# APPLICATION MODEL INFORMATION \n%s\n\n"%(separ,separ))     
        i.write("BEGIN MODEL_INPUT_FILES	KEYWORDS\n")
        i.write("  modinfile=%smodin.sen\n"%project)
        i.write("  templatefile=%stemplin.tpl\n"%project)
        i.write("END MODEL_INPUT_FILES\n\n")
        
        i.write("BEGIN MODEL_OUTPUT_FILES  KEYWORDS\n")
        i.write("  modoutfile=%smodout._os\n"%project)
        i.write("  instructionfile=%sobs.instructions\n"%project)
        i.write("  category=obs\n")
        i.write("END MODEL_OUTPUT_FILES")

        i.close()

        ## if flags g stop else otherwise run UCODE##    
        #if flags["g"]:
        if (mode=='config'):
            # Remove sampling points vector copy ##
            grass.run_command("g.remove",vect="%s"%points_copy) 
            
        else:
            ## Set UCODE command arguments ##
            filein = "%sucode.in"%project
            if len(prefix)>0:
                fileout = "%s/%s"%(project,prefix)
            else:
                fileout = "%sout"%project
                
            if len(exec_ucode)==0:
                exec_ucode="ucode"
                
            ## Run UCODE ##  
            proc00 = subprocess.Popen(["%s"%exec_ucode, "%s"%filein, "%s"%fileout],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            res00 = proc00.communicate()[0].split('\n')

            # Remove sampling points vector copy ##
            grass.run_command("g.remove",vect="%s"%points_copy) 


    #if flags['r']:
    if (mode=='advance'):
        ## search for .in file ##
        import os
        files_in=[]
        os.chdir("%s"%project)
        for files in os.listdir("."):
            if files.endswith(".in"):
                files_in.append(files)
                
        if len(files_in)==1:
            filein="%s%s"%(project,files_in[0])
        elif len(files_in)>1:            
           grass.fatal(_('More .in files found for the specified folder path'))
        else:
            grass.fatal(_('No .in file found for the specified folder path'))

        if len(prefix)>0:
            fileout = "%s/%s"%(project,prefix)
        else:
            fileout = "%sout"%project 
             
        if len(exec_ucode)==0:
            exec_ucode="ucode"      
                  
        ## Run UCODE ##                       
        proc00 = subprocess.Popen(["%s"%exec_ucode, "%sucode.in"%project, "%s"%fileout],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        res00 = proc00.communicate()[0].split('\n')
    
    
if __name__ == "__main__":
    options,flags = grass.parser()
    sys.exit(main())
