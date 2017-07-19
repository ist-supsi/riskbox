#!/usr/bin/env python
#
############################################################################
#
# MODULE:      trigrs.py
# AUTHOR(S):   Unknown; updated to GRASS 5.7 by Michael Barton
#              Converted to Python by Glynn Clements
# PURPOSE:     Uses d.his to drape a color raster over a shaded relief map
# COPYRIGHT:   (C) 2004,2008,2009 by the GRASS Development Team
# 
#              This program is free software under the GNU General Public
#              License (>=v2). Read the file COPYING that comes with GRASS
#              for details.
#
#############################################################################%Module
#%Module
#% description: Shallow landslide vulnerability assessment
#% keywords: raster
#%End
#%flag
#% key:r
#% description: set the region project
#% guisection: Required
#%End
#%option
#% key: project
#% type: string
#% description: project path
#% required: yes
#%end
#%option
#% key: mmax
#% type: integer
#% description: Maximum number of terms in series solution for saturated infiltration
#% answer: 20 
#% guisection: Required
#%end
#%option
#% key: nmax
#% type: integer
#% description: Maximum number of roots in series solution for unsaturated infiltration
#% answer: 30
#% guisection: Required
#%end
#%option
#% key: tx
#% type: integer
#% description: Time steps multiplier
#% guisection: Required
#%end
#%option
#% key: nzs
#% type: integer
#% description: Number of vertical increments
#% guisection: Required
#%end
#%option
#% key: min_slope
#% type: double
#% answer: 0.0
#% description: Minimum slope angle
#% guisection: Required
#%end
#%option
#% key: zmin
#% type: double
#% answer: 0.01
#% description: Shallowest depth at which to compute pore pressure and SF
#% guisection: Required
#%end
#%option
#% key: t
#% type: double
#% description: Computational time length [s]
#% guisection: Required
#%end
#%option
#% key: times
#% type:string
#% description: Times of output raster (t1,t2...)
#% guisection: Output
#%end
#%option
#% key: runoff
#% type:string
#% gisprompt: new,cell,raster
#% description: Runoff raster map
#% guisection: Output
#%end
#%option
#% key: minsf
#% type:string
#% gisprompt: new,cell,raster
#% description: Minimum safety factor raster map
#% guisection: Output
#%end
#%option
#% key: z_minsf
#% type:string
#% gisprompt: new,cell,raster
#% description: Depth of minimum safety factor raster map
#% guisection: Output
#%end
#%option
#% key: p_minsf
#% type:string
#% gisprompt: new,cell,raster
#% description: Pore pressure of minimum safety factor raster map
#% guisection: Output
#%end
#%option
#% key: infiltr_rate
#% type:string
#% gisprompt: new,cell,raster
#% description: Actual infiltration rate raster map
#% guisection: Output
#%end
#%option
#% key: basal_flux
#% type:string
#% gisprompt: new,cell,raster
#% description: Unsaturated zone basal flux raster map
#% guisection: Output
#%end
#%option
#% key: list_out
#% type: string
#% gisprompt: new_file,file,output
#% key_desc: name
#% description: Listing of pressure head and factor safety
#% guisection: Output
#%end
#%flag
#% key:a
#% description: use analytic solution for fillable porosity
#% guisection: Advanced
#%End
#%flag
#% key:m
#% description: log mass balance result
#% guisection: Advanced
#%End
#%option
#% key: flow_direction
#% type:string
#% description: Flow direction type
#% options: gener,slope,hydro
#% answer: gener
#% guisection: Advanced
#%end
#%option
#%  key: psi
#%  type: string
#%  options:the initial water table, the capillary fringe
#%  answer:the initial water table
#%  description: Use of the insaturated computations above
#% guisection: Advanced
#%end

import sys,os,re,tempfile,math,time,shutil,re,fileinput
space = re.compile(r'\s+')
multiSpace = re.compile(r"\s\s+")
import grass.script as grass
import subprocess
from subprocess import Popen

def main():
    proj = options['project']
    mmax = options['mmax']
    nmax = options['nmax']
    nzs = options['nzs']
    tx = options['tx']
    t = options['t']
    min_slope = options['min_slope']
    zmin = options['zmin']    
    runoff_map = options['runoff']
    minsf = options['minsf']
    z_minsf = options['z_minsf']    
    p_minsf = options['p_minsf']   
    infiltr_rate = options['infiltr_rate']   
    basal_flux = options['basal_flux']   
    times = options['times']    
    flow = options['flow_direction']    
    psi = options['psi']  
    list_out = options['list_out']  




######### tools ############

    def InfoCurrentRegion():
        region = grass.region()
        north = [region[key] for key in "nsew"][0]
        sud = [region[key] for key in "nsew"][1]
        est = [region[key] for key in "nsew"][2]
        west = [region[key] for key in "nsew"][3]
        nsres =  (region['nsres'])
        ewres =  (region['ewres'])
        cols = (region['cols'])
        rows = (region['rows'])
        return north,sud,west,est,nsres,ewres,rows,cols   
    
    def CheckVar(input_data,name):
        if input_data:
            pass
        else:
            grass.fatal('A %s value is required'%name)
    
    
    def SetOut(map):            
        if len(map)==0:
            check_map='F'
        else:
            check_map='T'
        return check_map
    
    
################## check input data #########################

    if len(proj)==0:
        grass.fatal('A project path is required')
    else:
        if os.path.exists(proj):
            pass
        else:
            grass.fatal('The project specified does not exist')        
            
            
    CheckVar(tx,'tx')        
    CheckVar(nzs,'nzs')
    CheckVar(t,'t')
    
    if len(times)==0:
        grass.fatal('Parameter times is required')
        
        
################# import projlog.log #######################

    print proj
    sys.path.append('%s'%proj)
    import projlog
################## create new_folder ######################
    if not os.path.exists(proj+'/TrigrsRes'):
        os.makedirs(proj+'/TrigrsRes')
################## check region ############################

    if not flags['r']:
        (n,s,w,e,nsr,ewr,nr,nc) = InfoCurrentRegion()   

        if n==projlog.log['north'] and s==projlog.log['sud'] and w==projlog.log['west'] and e==projlog.log['est'] and nsr==projlog.log['nsres'] and ewr==projlog.log['ewres'] and nr==projlog.log['rows'] and nc==projlog.log['columns']:
            print 'region ok'
        else:
            grass.fatal ('The current region is not the project region: set the flag r')

################## modify region ############################

    if flags['r']:
        grass.run_command("g.region",n=projlog.log['north'],s=projlog.log['sud'],w=projlog.log['west'],e=projlog.log['est'],nsres=projlog.log['nsres'],ewres=projlog.log['ewres'])
        #grass.run_command("g.region","p")

################## extract TopoIndex data ####################

    lines = open(proj+'/'+'TopoIndexLog.txt', 'r').readlines()
    line = lines[-5]
    a = space.split(line.strip())
    imax=a[0]
    row=a[1]
    col=a[2]
    nwf=a[3]

################# extract n_times data ###################

    list_times = list(times.split(','))  
    n_times =len(list_times)  

################## output ################################

    check_runoff_map = SetOut(runoff_map)
    check_minsf = SetOut(minsf)
    check_z_minsf = SetOut(z_minsf)
    check_p_minsf = SetOut(p_minsf)
    check_infiltr_rate = SetOut(infiltr_rate)
    check_basal_flux = SetOut(basal_flux)
    check_list_out = SetOut(list_out)
    
    if flags['a']:
        check_fill = 'T'
    else:
        check_fill = 'F'
        
    if psi=='the initial water table':
        check_psi='F'
    else:
        check_psi='T'
        
    if flags['m']:
        check_mass = 'T'
    else:
        check_mass = 'F'

        
    
    print check_runoff_map
    print check_minsf
    print check_z_minsf
    print check_p_minsf
    print check_infiltr_rate
    print check_basal_flux
    print check_list_out
      
################## write tr_in  ##########################

    f = open(proj+'/tr_in.txt', 'w')
    f.write('Name of project (up to 255 characters)\n')
    f.write('%s\n'%projlog.log['proj'])
    f.write('imax, row, col, nwf, tx, nmax\n')
    f.write('%s,\t%s,\t%s,\t%s,\t%s,\t%s\n'%(imax,row,col,nwf,tx,nmax))
    f.write('nzs, mmax, nper,  zmin,  uww,     t, zones\n')
    f.write('%s,\t%s,\t%s,\t%s,\t9.8e3,\t%s,\t%s\n'%(nzs,mmax,projlog.log['n_per'],zmin,t,projlog.log['n_zone']))
    f.write('zmax,   depth,   rizero,  Min_Slope_Angle (degrees)\n')
    f.write('-3.001,\t-2.4,\t-1.0e-9,\t%s\n'%min_slope)
    for x in projlog.log['range_cat_zone']:
        f.write('zone,\t%s\n'%x)
        f.write('cohesion,phi,  uws,   diffus,   K-sat, Theta-sat,Theta-res,Alpha\n')
        f.write(projlog.log['cohes_%s'%x]+',\t'+projlog.log['phi_%s'%x]+',\t'+projlog.log['uws_%s'%x]+',\t'+projlog.log['diffus_%s'%x]+',\t'+projlog.log['ksat_%s'%x]+',\t'+projlog.log['thetasat_%s'%x]+',\t'+projlog.log['thetares_%s'%x]+',\t'+projlog.log['a_%s'%x]+'\n')
    f.write('cri(1), cri(2), ..., cri(nper)\n')
    y=1
    while y<projlog.log['n_per']:
        f.write('-9.e-5,\t')
        y+=1
    f.write('-3.e-7\n')
    f.write('capt(1), capt(2), ..., capt(n), capt(n+1)\n')
    f.write('0,\t'+',\t'.join(projlog.log['list_capt'])+'\n')
    f.write('File name of slope angle grid (slofil)\n')
    f.write(projlog.log['slope_path']+'\n')
    f.write('File name of property zone grid (zonfil)\n')
    f.write(projlog.log['vzones_path']+'\n')
    f.write('File name of depth grid (zfil)\n')
    f.write(projlog.log['zmax_path']+'\n')
    f.write('File name of initial depth of water table grid   (depfil)\n')
    f.write(projlog.log['depthwt_path']+'\n')
    f.write('File name of initial infiltration rate grid   (rizerofil)\n')
    f.write(projlog.log['rizero_path']+'\n')
    f.write('List of file name(s) of rainfall intensity for each period, (rifil())\n')
    for item in projlog.log['rifil_path']:
        f.write(item+'\n')
    f.write('File name of grid of D8 runoff receptor cell numbers (nxtfil)\n')
    f.write(proj+'/'+'TIdscelGrid_%s.asc\n'%projlog.log['proj'])
    print proj+'/'+'TIdscelGrid_%s.asc\n'%projlog.log['proj']
    f.write('File name of list of defining runoff computation order (ndxfil)\n')
    f.write(proj+'/'+'TIcelindxList_%s.txt\n'%projlog.log['proj'])
    f.write('File name of list of all runoff receptor cells  (dscfil)\n')
    f.write(proj+'/'+'TIdscelList_%s.txt\n'%projlog.log['proj'])
    f.write('File name of list of runoff weighting factors  (wffil)\n')
    f.write(proj+'/'+'TIwfactorList_%s.txt\n'%projlog.log['proj'])
    f.write('Folder where output grid files will be stored  (folder)\n')
    f.write(proj+'/TrigrsRes/\n')
    f.write('Identification code to be added to names of output files (suffix)\n')
    f.write(projlog.log['proj']+'\n')
    f.write('Save grid files of runoff? Enter T (.true.) or F (.false.)\n')
    f.write(check_runoff_map+'\n')
    f.write('Save grid of minimum factor of safety? Enter T (.true.) or F (.false.)\n')
    f.write(check_minsf+'\n')
    f.write('Save grid of depth of minimum factor of safety? Enter T (.true.) or F (.false.)\n')
    f.write(check_z_minsf+'\n')
    f.write('Save grid of pore pressure at depth of minimum factor of safety? Enter T (.true.) or F (.false.)\n')
    f.write(check_p_minsf+'\n')
    f.write('Save grid files of actual infiltration rate? Enter T (.true.) or F (.false.)\n')
    f.write(check_infiltr_rate+'\n')
    f.write('Save grid files of unsaturated zone basal flux? Enter T (.true.) or F (.false.)\n')
    f.write(check_basal_flux+'\n')
    f.write('Save listing of pressure head and factor of safety ("flag")? (Enter -2 detailed, -1 normal, 0 none)\n')
    if check_list_out=='T':
        f.write('-2\n')
    else:
        f.write('0\n')
    f.write('Number of times to save output grids\n')
    f.write('%s\n'%n_times)
    f.write('Times of output grids\n')
    f.write(',\t'.join(list_times)+'\n')
    f.write('Skip other timesteps? Enter T (.true.) or F (.false.)\n')
    f.write('F\n')
    f.write('Use analytic solution for fillable porosity?  Enter T (.true.) or F (.false.)\n')
    f.write(check_fill+'\n')
    f.write('Estimate positive pressure head in rising water table zone (i.e. in lower part of unsat zone)?  Enter T (.true.) or F (.false.)\n')
    f.write('T\n')
    f.write('Use psi0=-1/alpha? Enter T (.true.) or F (.false.) (False selects the default value, psi0=0)\n')
    f.write(check_psi+'\n')
    f.write('Log mass balance results?   Enter T (.true.) or F (.false.)\n')
    f.write(check_mass+'\n')
    f.write('Flow direction (enter "gener", "slope", or "hydro")\n')
    f.write(flow+'\n')   
    f.close()


################## launch trigrs  ##########################


    folder0 = proj+'/'
    proc00 = subprocess.Popen(["/usr/local/Trigrs/Trigrs", "%s"%folder0],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    res00 = proc00.communicate()[0].split('\n')

################## import trigrs outputs ##########################
    respath=proj+'/TrigrsRes'
    if os.path.exists(respath):
        dirList=os.listdir(respath)
        filelist=[]
        for fname in dirList:
            if not fname.split("_")[0] == "TRnon":
                filelist.append(fname)
            else:
                grass.fatal("Error: solution non-convergent, consider adjust parameters")

        print filelist
        
        if len(filelist)<1:
            grass.fatal("Error: no output generated")
            
        filelist.sort()
               
        """
        if flags["-o"]:
            ov = True
        else:
            ov = False
        """
        
        for f in filelist:
            pref = f.split("_")
            if pref[0] == "TRfs":
                if pref[-1].split(".")[0].isdigit():
                    print int(pref[-1].split(".")[0])
                    print n_times
                    grass.run_command('r.in.arc',input=respath + "/" + f,output=minsf+"_"+list_times[int(pref[-1].split(".")[0])-1],overwrite = True)
                else:
                    grass.run_command('r.in.arc',input= respath + "/" + f,output=minsf,overwrite = True)
            if pref[0] == "TRz":
                if pref[-1].split(".")[0].isdigit():
                    grass.run_command('r.in.arc',input=respath + "/" + f,output=z_minsf+"_"+list_times[int(pref[-1].split(".")[0])-1],overwrite = True)
                else:
                    grass.run_command('r.in.arc',input=respath + "/" + f,output=z_minsf,overwrite = True)
            if pref[0] == "TRp":
                if pref[-1].split(".")[0].isdigit():
                    grass.run_command('r.in.arc',input=respath + "/" + f,output=p_minsf+"_"+list_times[int(pref[-1].split(".")[0])-1],overwrite = True)
                else:
                    grass.run_command('r.in.arc',input=respath + "/" + f,output=p_minsf,overwrite = True)



    
######################## remove trigrs result files ###########################################    

    for item in dirList:
        os.remove(respath+'/'+item)

#####################################################################


if __name__ == "__main__":
    options, flags = grass.parser()
    main()
