#!/usr/bin/env python
#
############################################################################
#
# MODULE:      trigrs.init.py
# AUTHOR(S):   Unknown; updated to GRASS 5.7 by Michael Barton
#              Converted to Python by Glynn Clements
# PURPOSE:     Uses d.his to drape a color raster over a shaded relief map
# COPYRIGHT:   (C) 2004,2008,2009 by the GRASS Development Team
# 
#              This program is free software under the GNU General Public
#              License (>=v2). Read the file COPYING that comes with GRASS
#              for details.
#
############################################################################
 
#%Module
#% description: Shallow landslide vulnerability assessment
#% keywords: raster
#%End
#%flag
#% key:p
#% description: print last parameters setting
#% guisection: Required
#%End
#%flag
#% key:m
#% description: modify an existent project
#% guisection: Required
#%End
#%option
#% key: project
#% type: string
#% description: Name of the project
#% required: yes
#%end
#%option
#% key: path
#% type: string
#% description: Folder path to create
#% required: yes
#% answer : /tmp
#%end
#%option
#% key: runoff 
#% type: integer
#% description: Weighting factor for runoff distribution
#% guisection: Runoff
#% options: -1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21
#%end
#%option
#% key: i
#% type: integer
#% description: Maximum number of iterations for cells indexing
#% guisection: Runoff
#%end
#%option
#% key: elevation
#% type: string
#% gisprompt: old,cell,raster
#% description: Elevation map
#% guisection: HydroGeoData
#%End
#%option
#% key: zmax
#% type: string
#% gisprompt: old,cell,raster
#% description: Basal boundary depth map
#% guisection: HydroGeoData
#%End
#%option
#% key: depthwt
#% type: string
#% gisprompt: old,cell,raster
#% description: Water table initial depth map
#% guisection: HydroGeoData
#%End
#%option
#% key: rizero
#% type: string
#% gisprompt: old,cell,raster
#% description: Steady background infiltration rate map
#% guisection: HydroGeoData
#%End
#%option
#% key: vzones
#% type: string
#% gisprompt: old,vector,vector
#% description: Hydrogeological zones map 
#% guisection: HydroGeoData
#%End
#%option
#% key: rifil
#% type: string
#% gisprompt: old,cell,raster
#% multiple: yes
#% description: rainfall intensity raster maps (ri1,ri2...)
#% guisection: MeteoData
#%End
#%option
#% key: capt
#% type: string
#% description: cumulative duration (in seconds) of rainfall intensity (capt1,capt2...)
#% guisection: MeteoData
#%End

import sys,os,re,tempfile,math,time,shutil,re,fileinput
import grass.script as grass
import subprocess
from subprocess import Popen

def main():
    proj = options['project']
    path = options['path']
    dtm = options['elevation']
    wf = options['runoff']
    it = options['i']
    zmax = options['zmax']   
    depthwt = options['depthwt']
    rizero = options['rizero']
    vzones = options['vzones']
    rifil = options['rifil']
    capt = options['capt']

    p=re.compile('\W')
    


    def SetMapsetMap(map):
        if len(map)!=0:
            c=map.find('@')
            if c==-1:
                env = grass.gisenv()
                map = map+'@'+env['MAPSET']
            else:
                map=map
        return map             
    
    dtm=SetMapsetMap(dtm)
    zmax=SetMapsetMap(zmax)
    depthwt=SetMapsetMap(depthwt)
    rizero=SetMapsetMap(rizero)
    vzones=SetMapsetMap(vzones)

    list_rifil= list(rifil.split(','))
    j=0
    for item in list_rifil:
        item_n=SetMapsetMap(item)
        list_rifil[j]=item_n
        j+=1
    rifil = ','.join(list_rifil)

################## tools ########################### 
       
    
    def Dict(key,value):
        log[key]=value    
    
    def ModDict(key,value):
        projlog.log[key]=value    
        
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
    
    def GetCtime(map):
        map_ok = map.split('@')[0]
        env = grass.gisenv()
        ctime = os.stat (env['GISDBASE']+env['LOCATION_NAME']+'/'+env['MAPSET']+'/'+'cellhd/'+map_ok).st_ctime
        return ctime
    
    def GetCtimeVect(map):
        map_ok = map.split('@')[0]
        env = grass.gisenv()
        ctime = os.stat (env['GISDBASE']+env['LOCATION_NAME']+'/'+env['MAPSET']+'/'+'vector/'+map_ok+'/head').st_ctime
        return ctime
        
    def GetDirection (input_map): 
        grass.run_command ("r.fill.dir", input=input_map, elevation=t_outfill, direction=t_direction, overwrite=True,quiet=True)
        grass.run_command("r.mapcalc", expression='%s = if (%s==135,1,if(%s==90,2,if(%s==45,3,if(%s==180,4,if(%s==360,6,if(%s==225,7,if(%s==270,8,if(%s==315,9,null()))))))))'%(direction,t_direction,t_direction,t_direction,t_direction,t_direction,t_direction,t_direction,t_direction),quiet=True) 

        
    def GetSlope(input_map):
        grass.run_command("r.slope.aspect", elevation=input_map, slope=slope, overwrite=True, quiet=True)     
        
    
    def GetIterations(input_i,row,col):
        if (input_i=='-1'):
            n_iter = int(math.ceil((1+((float(col)+float(row))/20))))
        else:
            n_iter = input_i
        return n_iter   
    
    def SetRegion(nn,ss,ww,ee):
        grass.run_command("g.region",n=nn,s=ss,w=ww,e=ee)
        
    def ExtractData(input_map,folder_map):
        grass.run_command("r.out.arc",input=input_map,output=folder_map,quiet=True)
        
    def RemDtmTemp ():
        grass.run_command("g.remove",rast=slope, quiet=True)
        grass.run_command("g.remove",rast=t_direction, quiet=True)
        grass.run_command("g.remove",rast=direction, quiet=True)
        grass.run_command("g.remove",rast=t_outfill, quiet=True)
        
    def Tpx_in():
        f = open(log['folder_proj']+'/tpx_in.txt', 'w')
        f.write('Name of project (up to 255 characters)\n')
        f.write('%s\n'%log['proj'])
        f.write('Rows, Columns, flow-direction numbering scheme (ESRI=1, TopoIndex=2)\n')
        f.write('%s,\t%s,\t2\n'%(log['rows'],log['columns']))
        f.write('Exponent, Number of iterations\n')
        f.write('%s,\t%s\n'%(log['runoff'],log['iterations']))
        f.write('Name of elevation grid file\n')
        f.write('%s\n'%log['dtm_path'])
        f.write('Name of direction grid\n')
        f.write('%s\n'%log['direction_path'])
        f.write('Save listing of D8 downslope neighbor cells?  Enter T (.true.) or F (.false.)\n')
        f.write('F\n')
        f.write('Save grid of D8 downslope neighbor cells? Enter T (.true.) or F (.false.)\n')
        f.write('T\n')
        f.write('Save cell index number grid?  Enter T (.true.) or F (.false.)\n')
        f.write('F\n')
        f.write('Save list of cell number and corresponding index number? Enter T (.true.) or F (.false.)\n')
        f.write('T\n')
        f.write('Save flow-direction grid remapped from ESRI to TopoIndex? Enter T (.true.) or F (.false.)\n')
        f.write('F\n')
        f.write('Name of folder to store output?\n')
        f.write('%s/\n'%log['folder_proj'])
        f.write('ID code for output files? (8 characters or less)\n')
        f.write('%s\n'%log['proj'])
        f.close()
        
    def Tpx_in_m():
        f = open(projlog.log['folder_proj']+'/tpx_in.txt', 'w')
        f.write('Name of project (up to 255 characters)\n')
        f.write('%s\n'%projlog.log['proj'])
        f.write('Rows, Columns, flow-direction numbering scheme (ESRI=1, TopoIndex=2)\n')
        f.write('%s,\t%s,\t2\n'%(projlog.log['rows'],projlog.log['columns']))
        f.write('Exponent, Number of iterations\n')
        f.write('%s,\t%s\n'%(projlog.log['runoff'],projlog.log['iterations']))
        f.write('Name of elevation grid file\n')
        f.write('%s\n'%projlog.log['dtm_path'])
        f.write('Name of direction grid\n')
        f.write('%s\n'%projlog.log['direction_path'])
        f.write('Save listing of D8 downslope neighbor cells?  Enter T (.true.) or F (.false.)\n')
        f.write('F\n')
        f.write('Save grid of D8 downslope neighbor cells? Enter T (.true.) or F (.false.)\n')
        f.write('T\n')
        f.write('Save cell index number grid?  Enter T (.true.) or F (.false.)\n')
        f.write('F\n')
        f.write('Save list of cell number and corresponding index number? Enter T (.true.) or F (.false.)\n')
        f.write('T\n')
        f.write('Save flow-direction grid remapped from ESRI to TopoIndex? Enter T (.true.) or F (.false.)\n')
        f.write('F\n')
        f.write('Name of folder to store output?\n')
        f.write('%s/\n'%projlog.log['folder_proj'])
        f.write('ID code for output files? (8 characters or less)\n')
        f.write('%s\n'%projlog.log['proj'])
        f.close()    
         
    def Tpx_Launch():
        folder0 = log['folder_proj']+'/'
        proc00 = subprocess.Popen(["/usr/local/Trigrs/TopoIndex", "%s"%folder0],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        res00 = proc00.communicate()[0].split('\n')
        last_line = res00[len(res00)-2]
        return last_line      
    
    def Tpx_Launch_m():
        folder0 = projlog.log['folder_proj']+'/'
        proc00 = subprocess.Popen(["/usr/local/Trigrs/TopoIndex", "%s"%folder0],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        res00 = proc00.communicate()[0].split('\n')
        last_line = res00[len(res00)-2]
        return last_line          
    

    def InitCheckMap(map_list,input_data,mess):
        if input_data in map_list:
            return
        else:
            grass.fatal('To create a new project %s is needed!'%mess)

            
    def CheckVar(name):                          
        if not name:
            check_name = 'False'
        else:
            check_name = 'True'
        return check_name


    def SetCheckValue(new_map,old_map,crt,maptype):
        if len(new_map)==0:
            check_map = 'False'
        else:
            if new_map==old_map:
                if maptype=='rast':
                    new_crt = GetCtime(new_map)
                if maptype=='vect':
                    new_crt = GetCtimeVect(new_map)
                if crt==new_crt:
                    check_map='False'
                else:
                    check_map='True'
            else:
                # verificare se esiste
                if maptype=='rast':
                    ex_rast = grass.list_strings('rast')
                    if new_map in ex_rast:
                        check_map='True'
                    else:
                        grass.fatal('%s does not exist'%new_map)
                if maptype=='vect':
                    ex_vect = grass.list_strings('vect')
                    if new_map in ex_vect:
                        check_map='True'
                    else:
                        grass.fatal('%s does not exist'%new_map)               
        return check_map
    
                
    def CheckCaptType(list_capt):
        for item in list_capt:
            try:
                float(item)
            except:
                grass.fatal('capt must be a numeric value')
                
                
    def CheckCaptValue(list_capt):  
        if len(list_capt)==1:
            pass
        else:
            for x in range(len(list_capt)-1):
                if float(list_capt[x]) < float(list_capt[x+1]):
                    pass
                else:
                    grass.fatal('capt(n) must be smaller than capt(n+1)')
            
                        
####################### Only print data setting ######################


    if flags['p']:
        sys.path.append(path+'/'+proj)
        import projlog
        print projlog.log

        print '##### LAST PARAMETERS #####'
        print '--------- Project ----------'
        print 'project_folder = %s'%(projlog.log['folder_proj'])
        print 'project name = %s'%(projlog.log['proj'])
        print '--------- Region ----------'
        print 'north = %s' %(projlog.log['north'])
        print 'south = %s' %(projlog.log['sud'])
        print 'east = %s' %(projlog.log['est'])
        print 'west = %s' %(projlog.log['west'])
        print 'nsres = %s' %(projlog.log['nsres'])
        print 'ewres = %s' %(projlog.log['ewres'])
        print 'rows = %s' %(projlog.log['rows'])
        print 'columns = %s' %(projlog.log['columns'])
        print '--------- Tpx data  ----------'
        print 'runoff = %s' %(projlog.log['runoff'])
        print 'it = %s' %(projlog.log['it'])
        print 'iterations = %s' %(projlog.log['iterations'])
        print '--------- Dtm map ----------'
        print 'name = %s' %(projlog.log['dtm'])
        print 'path = %s' %(projlog.log['dtm_path'])
        print 'ctime = %s' %(projlog.log['ctime_dtm'])
        print '--------- Direction map ----------'
        print 'path = %s' %(projlog.log['direction_path'])
        print '--------- Slope map ----------'
        print 'path = %s' %(projlog.log['slope_path'])
        print '--------- Basal boundary depth map ----------'
        print 'name = %s' %(projlog.log['zmax'])
        print 'path = %s' %(projlog.log['zmax_path'])
        print 'ctime = %s' %(projlog.log['ctime_zmax'])
        print '--------- Water table initial depth map ----------'
        print 'name = %s' %(projlog.log['depthwt'])
        print 'path = %s' %(projlog.log['depthwt_path'])
        print 'ctime = %s' %(projlog.log['ctime_depthwt'])
        print '--------- Steady background infiltration rate map ----------'
        print 'name = %s' %(projlog.log['rizero'])
        print 'path = %s' %(projlog.log['rizero_path'])
        print 'ctime = %s' %(projlog.log['ctime_rizero'])
        print '--------- Hydrogeological data ----------'
        print 'vector_name=%s'%(projlog.log['vzones'])
        print 'vector_ctime=%s'%(projlog.log['ctime_vzones'])
        print 'raster_path=%s'%(projlog.log['vzones_path'])
        print 'n_zone=%s'%(projlog.log['n_zone'])
        for x in projlog.log['range_cat_zone']:
            print '--------zone%s--------'%x
            print 'cohesion='+projlog.log['cohes_%s'%x]
            print 'phi='+projlog.log['phi_%s'%x]
            print 'uws='+projlog.log['uws_%s'%x]
            print 'diffus='+projlog.log['diffus_%s'%x]
            print 'ksat='+projlog.log['ksat_%s'%x]
            print 'thetasat='+projlog.log['thetasat_%s'%x]
            print 'thetares='+projlog.log['thetares_%s'%x]
            print 'a='+projlog.log['a_%s'%x]
        print '--------- Meteo data ----------'
        print 'n_periods=%s'%projlog.log['n_per']
        for x in range(projlog.log['n_per']):
            print '--------period_%s--------'%(x+1)
            print 'name=%s'%projlog.log['list_rifil'][x]
            print 'ctime=%s'%projlog.log['list_ctime_rifil'][x]  
            print 'path=%s'%projlog.log['rifil_path'] [x]
            print 'cumul_duration=%s'%projlog.log['list_capt'][x]
            
        sys.exit()


########################## Data check #################################

    if not flags['m']:        
        ex_rast = grass.list_strings('rast')
        ex_vect = grass.list_strings('vect')
        
        check_wf = CheckVar(wf)
        if check_wf=='False':
            grass.fatal('To create a new project runoff weighting factor is needed')      
            
        check_it = CheckVar(it)
        if check_it=='False':
            grass.fatal('To create a new project iteration value is needed')  

        InitCheckMap(ex_rast,dtm,'elevation map')
        InitCheckMap(ex_rast,zmax,'basal boundary depth map')
        InitCheckMap(ex_rast,depthwt,'water table initial depth map')
        InitCheckMap(ex_rast,rizero,'steady background infiltration rate map')
        InitCheckMap(ex_vect,vzones,'hydrogeological zones map')
      
       ################################# meteo data ######################################
                
        if len(rifil)==0:
            grass.fatal('To create a new project at least one map of intensity rainfall is required')
        else:
            list_rifil= list(rifil.split(','))
            for item in list_rifil:
                if item in ex_rast:
                    pass
                else:
                    grass.fatal('%s does not exist'%item)
            
        if len(capt)==0:
            grass.fatal('To create a new project cumulative duration of each rainfall intensity is required')   
        else:
            list_capt= list(capt.split(','))


            # Check vari #
            if len(list_capt)!=len(list_rifil):                                      
                grass.fatal('intensity and duration data does not correspond')
              
            CheckCaptType(list_capt)
            CheckCaptValue(list_capt)



################ create/modify Project e init dictionary ######################            

    if not flags['m']:
        print 'create project'
        if os.path.exists(path):
            folder = os.path.normpath(path+'/'+proj)
            if os.access(path,os.W_OK):
                os.makedirs(folder)
                log={}
                Dict('folder_proj',folder)
                Dict('proj',proj)
            else:
                grass.fatal ('no permissions to write in the specified path.')
        else:
            grass.fatal ('the specified path does not exist.')   

    else:
        if os.path.exists(path+'/'+proj):
            print 'modify project'
            pass
        else:
            grass.fatal('The specified project does not exist.')            
            
############################# Verify input value modified ############################

    if flags['m']:        
        sys.path.append(path+'/'+proj)
        import projlog
        check_dtm = SetCheckValue(dtm,projlog.log['dtm'],projlog.log['ctime_dtm'],'rast')
        check_zmax = SetCheckValue(zmax,projlog.log['zmax'],projlog.log['ctime_zmax'],'rast')
        check_depthwt = SetCheckValue(depthwt,projlog.log['depthwt'],projlog.log['ctime_depthwt'],'rast')
        check_rizero = SetCheckValue(rizero,projlog.log['rizero'],projlog.log['ctime_rizero'],'rast')
        #check_vzones = SetCheckValue(vzones,projlog.log['vzones'],projlog.log['ctime_vzones'],'vect')
   
        ############# vzones provvisorio ##############
        
        if len(vzones)==0:
            check_vzones='False'
        else:
            ex_vect = grass.list_strings('vect')
            if vzones in ex_vect:
                check_vzones='True'
            else:
                grass.fatal('%s map does not exist'%vzones)
        
   
       ##################################################
        
        check_wf=CheckVar(wf)
        if check_wf=='True' and wf==projlog.log['runoff']:
            check_wf='False'
        else:
            pass
        
        check_it=CheckVar(it)
        if check_it=='True' and it==projlog.log['it']:
            check_it='False'
        else:      
            pass  

        ############ check meteo data #############
        
        
        if len(rifil)==0:
            check_ri = 'False'
        else:
            new_list_rifil= list(rifil.split(','))
            if new_list_rifil==projlog.log['list_rifil']:
                nct=[]
                for item in new_list_rifil:
                    new_ct = GetCtime(item)
                    nct.append(new_ct)
                if nct==projlog.log['list_ctime_rifil']:
                    check_ri = 'False'
                else:
                    check_ri = 'True'
            else:
                ex_rast = grass.list_strings('rast')
                for item in new_list_rifil:
                    if item in ex_rast:
                        check_ri = 'True'
                    else:
                        grass.fatal('%s does not exist'%item)
                    
        if len(capt)==0:
            check_capt='False'
        else:
            new_list_capt= list(capt.split(','))
            print new_list_capt
            if len(new_list_capt)!=len(new_list_rifil): 
                grass.fatal('intensity and duration data does not correspond')
            CheckCaptType(new_list_capt)
            CheckCaptValue(new_list_capt)
            
            if new_list_capt == projlog.log['list_capt']:
                check_capt = 'False'
            else:     
                check_capt='True'
                    

        print 'dtm',check_dtm       
        print 'zmax',check_zmax
        print 'depthwt',check_depthwt
        print 'rizero',check_rizero
        print 'wf',check_wf
        print 'it',check_it
        print 'vzones', check_vzones
        print 'rifil', check_ri
        print 'capt', check_capt
        

########################## Insert input data in dictionary ################

    if not flags['m']:
        Dict("dtm",dtm)
        Dict("runoff",wf)
        Dict("it",it)
        Dict("zmax",zmax)
        Dict("depthwt",depthwt)
        Dict("rizero",rizero)
        Dict("vzones",vzones)
        
        ct_dtm = GetCtime(dtm)
        ct_zmax =  GetCtime(zmax)
        ct_depthwt = GetCtime(depthwt)
        ct_rizero = GetCtime(rizero)
        ct_vzones = GetCtimeVect(vzones)
        
        Dict("ctime_dtm",ct_dtm)
        Dict("ctime_zmax",ct_zmax)
        Dict("ctime_depthwt",ct_depthwt)
        Dict("ctime_rizero",ct_rizero)
        Dict("ctime_vzones",ct_vzones)
        
        
        # rifil #
        list_rifil= list(rifil.split(','))
        Dict('n_per',len(list_rifil))
        Dict('list_rifil',list_rifil)
        ct=[]
        for item in list_rifil:
            ct_ri = GetCtime(item)
            ct.append(ct_ri)
        Dict('list_ctime_rifil',ct)
        
            
        # capt # 
        list_capt= list(capt.split(','))
        Dict('list_capt',list_capt)

########################## Current region data ############################   

    if not flags['m']:     
        print 'extract current region value'   
        (n,s,w,e,nsr,ewr,nr,nc) = InfoCurrentRegion()   
        Dict('north', n)
        Dict('sud', s)
        Dict('west', w)
        Dict('est', e)
        Dict('nsres', nsr)
        Dict('ewres', ewr)
        Dict('rows', int(nr))
        Dict('columns', int(nc))

    if flags['m']:
        print 'verify current region'
        (n,s,w,e,nsr,ewr,nr,nc) = InfoCurrentRegion()  
        print n,s,w,e,nsr,ewr,nr,nc 
        
        if n==projlog.log['north'] and s==projlog.log['sud'] and w==projlog.log['west'] and e==projlog.log['est'] and nsr==projlog.log['nsres'] and ewr==projlog.log['ewres'] and nr==projlog.log['rows'] and nc==projlog.log['columns']:
            print 'region same'
        else:
            grass.fatal ('The current region is not the project region: modify or create a new project')
            
###############  Set DTM, slope and direction ###################


    if not flags['m']:
        print 'create dtm...'
        #  setup new temporary region  #
        new_n = n+nsr
        new_s = s-nsr
        new_w = w-ewr
        new_e = e+ewr
        SetRegion(new_n,new_s,new_w,new_e)
        
        #  create slope and direction  #   
        t_direction = p.sub("",tempfile.mktemp())
        t_outfill = p.sub("",tempfile.mktemp())
        direction = p.sub("",tempfile.mktemp())
        slope = p.sub("",tempfile.mktemp())
        GetDirection(dtm)
        GetSlope(dtm)
                
        #  set originally user region  #
        SetRegion(n,s,w,e)
                
        # log #
        env = grass.gisenv()
        Dict('dtm_path', log['folder_proj']+'/'+ str(dtm)+'.asc')
        Dict('direction_path', (log['folder_proj']+'/'+'direction@%s.asc')%env['MAPSET'])
        Dict('slope_path', (log['folder_proj']+'/'+'slope@%s.asc')%env['MAPSET'])
 
          
        #  extract data  #
        ExtractData(dtm,log['dtm_path'])
        ExtractData(direction,log['direction_path'])
        ExtractData(slope,log['slope_path'])    
        
        # remove temporary data #
        RemDtmTemp()    
    
    
    if flags['m'] and check_dtm=='False':
        print 'old dtm,slope,aspect'
        pass 
    
    
    
    if flags['m'] and check_dtm=='True':
        print 'cambiare dem'
        
        # delete old file #
        os.remove(projlog.log['dtm_path'])
        os.remove(projlog.log['direction_path'])
        os.remove(projlog.log['slope_path'])
        print 'dtm,slope,direction delete'
          
        #  setup new temporary region  #
        new_n = n+nsr
        new_s = s-nsr
        new_w = w-ewr
        new_e = e+ewr
        SetRegion(new_n,new_s,new_w,new_e)
        print 'new temporary region'
        
        #  create slope and direction  #   
        t_direction = p.sub("",tempfile.mktemp())
        t_outfill = p.sub("",tempfile.mktemp())
        direction = p.sub("",tempfile.mktemp())
        slope = p.sub("",tempfile.mktemp())
        GetDirection(dtm)
        GetSlope(dtm)  
        print 'create slope and direction'
        
        #  set originally user region  #
        SetRegion(n,s,w,e)
        print 'set originally user region'
        
        # modify dict #
        ModDict('dtm',dtm)       
        ct_dtm = GetCtime(dtm)
        ModDict('ctime_dtm',ct_dtm)   
        ModDict('dtm_path', projlog.log['folder_proj']+'/'+ str(dtm)+'.asc')     
        print projlog.log['dtm']
        print projlog.log['ctime_dtm']
        print projlog.log['dtm_path']
        
        # extract data #
        ExtractData(dtm,projlog.log['dtm_path'])
        ExtractData(direction,projlog.log['direction_path'])
        ExtractData(slope,projlog.log['slope_path'])  

  
        # remove temporary data #
        RemDtmTemp()                   


        
########################  Iterations  #############################

    if not flags['m']:
        print 'get iteration'
        iterations = GetIterations(it,nr,nc)
        print iterations
        Dict('iterations',iterations)  
        
    if flags['m'] and check_it=='False':
        print 'old it,iterations'
        pass
    
    if flags['m'] and check_it=='True':
        print 'new it,iterations'
        ModDict('it',it)
        iterations = GetIterations(it,nr,nc)
        ModDict('iterations',iterations)
        
       
##############################  Runoff ############################

    if not flags['m']:
        pass
    
    if flags['m'] and check_wf=='False':
        print 'old wf'
        pass
    
    if flags['m'] and check_wf=='True':
        print 'new wf'
        ModDict('runoff',wf)
    

############################# Tpx #################################

    if not flags['m']:
        print 'new_tpx'
        Tpx_in()        
        result=Tpx_Launch ()
        print result

        if (result == ' Correcting cell index numbers') and (it=='-1'):
            grass.fatal ('the index solution does not converge: a larger number of iterations is required\n.Default number of iterations is set to %s.\n If this error persist DEM is probably not hydrologically correct.'%iterations)
        if (result == ' Correcting cell index numbers')and (it!='-1'):
            grass.fatal ('the index solution does not converge: a larger number of iterations is required.\nIf this error persist DEM is probably not hydrologically correct.')    
    
    if flags['m'] and (check_wf=='True' or check_it=='True' or check_dtm=='True'):
        print 'rifo_tpx'

        Tpx_in_m()

        result=Tpx_Launch_m()
     
        if (result == ' Correcting cell index numbers') and (it==-1):
            grass.fatal ('the index solution does not converge: a larger number of iterations is required\n.Default number of iterations is set to %s.\n If this error persist DEM is probably not hydrologically correct.'%iterations)
        if (result == ' Correcting cell index numbers')and (it!=-1):
            grass.fatal ('the index solution does not converge: a larger number of iterations is required.\nIf this error persist DEM is probably not hydrologically correct.')            
      
      
    if flags['m'] and (check_wf=='False' and check_it=='False' and check_dtm=='False'):
        print 'vecchio_tpx'
        pass
    

         
######################## zmax map #################################

    if not flags['m']:    
        print 'create zmax.asc'
        Dict('zmax_path', log['folder_proj']+'/'+ str(zmax)+'.asc')
        ExtractData(zmax,log['zmax_path'])
        
    if flags['m'] and check_zmax=='False':
        print 'old zmax'
        pass   
        

    if flags['m'] and check_zmax=='True':
        
        # delete old file #
        os.remove(projlog.log['zmax_path'])
        print 'old_zmax removed'
        
        # modify Dict #
        ModDict('zmax',zmax)    
        ct_zmax = GetCtime(zmax)
        ModDict('ctime_zmax',ct_zmax)  
        ModDict('zmax_path', projlog.log['folder_proj']+'/'+ str(zmax)+'.asc')
 
        # extract new zmax #
        print 'new_z_max'
        ExtractData(zmax,projlog.log['zmax_path'])
        
            
#######################   depthwt map   ###########################

    if not flags['m']:    
        print 'create depthwt.asc'
        Dict('depthwt_path', log['folder_proj']+'/'+ str(depthwt)+'.asc')
        ExtractData(depthwt,log['depthwt_path'])
        
    if flags['m'] and check_depthwt=='False':
        print 'old depthwt'
        pass       
        
    if flags['m'] and check_depthwt=='True':
        
        # delete old file #
        os.remove(projlog.log['depthwt_path'])
        print 'old_depthwt removed'
        
        # modify Dict #
        ModDict('depthwt',depthwt)    
        ct_depthwt = GetCtime(depthwt)
        ModDict('ctime_depthwt',ct_depthwt)  
        ModDict('depthwt_path', projlog.log['folder_proj']+'/'+ str(depthwt)+'.asc')
 
        # extract new zmax #
        print 'new_depthwt'
        ExtractData(depthwt,projlog.log['depthwt_path'])      
        
          
#########################  rizero map  ################################  

    if not flags['m']:    
        print 'create rizero.asc'
        Dict('rizero_path', log['folder_proj']+'/'+ str(rizero)+'.asc')
        ExtractData(rizero,log['rizero_path'])

    if flags['m'] and check_rizero=='False':
        print 'old rizero'
        pass    
    
    if flags['m'] and check_rizero=='True':
        
        # delete old file #
        os.remove(projlog.log['rizero_path'])
        print 'old_rizero removed'
        
        # modify Dict #
        ModDict('rizero',rizero)    
        ct_rizero = GetCtime(rizero)
        ModDict('ctime_rizero',ct_rizero)  
        ModDict('rizero_path', projlog.log['folder_proj']+'/'+ str(rizero)+'.asc')
 
        # extract new rizero #
        print 'new_rizero'
        ExtractData(rizero,projlog.log['rizero_path'])      
    
        

#########################  zones map  ################################   
 
    if not flags['m']:
        # create raster zone #
        t_rast_zones = p.sub("",tempfile.mktemp())
        grass.run_command ("v.to.rast", input=vzones, output=t_rast_zones, use='cat', type='area', overwrite=True, quiet=True)
        
        # insert dict #
        env = grass.gisenv()       
        Dict('vzones_path', (log['folder_proj']+'/'+'%s.asc')%(vzones))
              
        # extract data #
        ExtractData(t_rast_zones,log['vzones_path'])    
      
        # remove temporary data #
        grass.run_command("g.remove",rast=t_rast_zones, quiet=True)
            
        # select cat and insert in dict#   
        res = grass.read_command("v.db.select",'c',map=vzones,columns='cat')
        range_cat = res.split('\n')
        range_cat.pop()
        Dict('range_cat_zone',range_cat)
        Dict('n_zone',len(range_cat))
        range_hdata = ['cohes','phi','uws','diffus','ksat','thetasat','thetares','a']
        for x in range_cat:
            for j in range_hdata:
                res = grass.read_command("v.db.select",'c',map=vzones,columns=j, where='cat=%s'%x)
                result = res.split('\n')[0]
                Dict('%s_%s'%(j,x),result)

    if flags['m'] and check_vzones=='False':
        print 'old vzones'
        pass         
        
    if flags['m'] and check_vzones=='True':
        print 'new_vzones'

        # delete old file #
        os.remove(projlog.log['vzones_path'])
        
        # delete old item dict #
        for x in projlog.log['range_cat_zone']:
            del projlog.log['cohes_%s'%x]
            del projlog.log['phi_%s'%x]
            del projlog.log['uws_%s'%x]
            del projlog.log['diffus_%s'%x]
            del projlog.log['ksat_%s'%x]
            del projlog.log['thetasat_%s'%x]
            del projlog.log['thetares_%s'%x]
            del projlog.log['a_%s'%x]
        print 'all removed'
         
        # create raster zone #
        t_rast_zones = p.sub("",tempfile.mktemp())
        grass.run_command ("v.to.rast", input=vzones, output=t_rast_zones, use='cat', type='area', overwrite=True, quiet=True)
        print 'new_raster_zone'
        
        # modify dict #
        ModDict('vzones',vzones) 
        ct_vzones = GetCtimeVect(vzones)
        ModDict('ctime_vzones',ct_vzones)  
        ModDict('vzones_path', projlog.log['folder_proj']+'/'+ str(vzones)+'.asc')
        
        # extract data #
        ExtractData(t_rast_zones,projlog.log['vzones_path'])        

        # remove temporary data #
        grass.run_command("g.remove",rast=t_rast_zones, quiet=True)              
        
        # select cat and insert new items dict #   
        res = grass.read_command("v.db.select",'c',map=vzones,columns='cat')
        range_cat = res.split('\n')
        range_cat.pop()                
        ModDict('range_cat_zone',range_cat)
        ModDict('n_zone',len(range_cat))

         
        range_hdata = ['cohes','phi','uws','diffus','ksat','thetasat','thetares','a']
        for x in range_cat:
            for j in range_hdata:
                res = grass.read_command("v.db.select",'c',map=vzones,columns=j, where='cat=%s'%x)
                result = res.split('\n')[0]
                projlog.log['%s_%s'%(j,x)]=result
       

#########################  rifil  ################################   

    if not flags['m']:
        print 'create rifil'
        list_rifil= list(rifil.split(','))
        path=[]
        for item in list_rifil:
            path_rifil=log['folder_proj']+'/'+ str(item)+'.asc'
            ExtractData(item,path_rifil) 
            path.append(path_rifil)           
        Dict('rifil_path', path)
        
        
    if flags['m'] and check_ri=='False':
        print 'old rifil'
        pass    
    
    
    if flags['m'] and check_ri=='True':
        print 'new rifil'
        # delete old file #
        for item in projlog.log['rifil_path']:
            os.remove(item)
        print 'all delete'
        
           
        # insert new items in dict and extract new_data #
        new_list_rifil= list(rifil.split(','))

        ModDict('list_rifil',new_list_rifil)
        ModDict('n_per',len(new_list_rifil))

        nct=[]
        npath = []
        for item in new_list_rifil:
            nct_ri = GetCtime(item)
            nct.append(nct_ri) 
            path_r=projlog.log['folder_proj']+'/'+ str(item)+'.asc'
            npath.append(path_r)
            ExtractData(item,path_r) 
        ModDict('list_ctime_rifil',nct)
        ModDict('rifil_path',npath)
        
#########################  capt  ##################################

    if not flags['m']:
        pass
    
    if flags['m'] and check_capt=='False':
        print 'old capt'
        pass  
    
    if flags['m'] and check_capt=='True':
        print 'new_capt'
        
        new_list_capt= list(capt.split(','))
        ModDict('list_capt',new_list_capt)        
        
#######################  Export dictionary  #######################

    if not flags['m']:
        file = open(log['folder_proj']+'/'+'projlog.py',"w")
        file.write('log=' + str (log))
        file.close()
        
        
    ## da aggiungere 'solo se ci sono modifiche'   
    if flags['m']:
        print 'new_dict'
        file = open(projlog.log['folder_proj']+'/'+'projlog.py',"w")
        file.write('log=' + str (projlog.log))
        file.close()

    
    
#####################################################################
    
    
if __name__ == "__main__":
    options, flags = grass.parser()
    main()
