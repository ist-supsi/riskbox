#!/usr/bin/env python
#

import sys
import os
import time
import grass.script as grass

#[value,variance,step]
s_rho    = [2100,2600,100]
s_vol    = [11000,14000,500]
s_width  = [20,30,1]
s_thick  = [20,30,1]
i_vel    = [30,45,2]
i_slope  = [30,45,2]
i_azimut = [200,250,5] 
i_depth  = [20,70,10]

"""
#python /home/ist/Desktop/r.impact.tsunami/r.impact.tsunami.py 
#elevation=dtmpatch lake=lake_depth s_rho=2500 s_vol=13000 s_width=26 
#s_thick=25 i_east=1527471.21021 i_north=5082144.5206 i_vel=35 i_slope=35 
#i_azimut=225 output=test_tsu i_depth=50 -ac --o
s_rho    = [2500,2501,100]
s_vol    = [13000,13001,500]
s_width  = [26,27,1]
s_thick  = [25,26,1]
i_vel    = [35,36,2]
i_slope  = [35,36,2]
i_azimut = [225,226,5] 
i_depth  = [60,70,10]
"""
out_prefix = "sens_"
out_report = "/home/ist/sensitivity.csv"

def getReport(mapname):
    ret = grass.read_command("r.univar",map=mapname,flags="g")
    info = {}
    for l in ret.split("\n"):
        s = l.split("=")
        if len(s)==2:
            info[s[0]]=float(s[1])
    return info

def main():
    index = ["loop"]
    pars=["rho","vol","width","thick","vel","slope","azimut","depth"]
    stats=["n","null_cells","cells","min","max","range","mean","mean_of_abs","stddev","variance","coeff_var","sum"]
    waves=["Hm","Tm","Lm","a","c","rm"]
    report=[]
    loop=0
    for rho in range(s_rho[0],s_rho[1],s_rho[2]):
        for vol in range(s_vol[0],s_vol[1],s_vol[2]):
            for width in range(s_width[0],s_width[1],s_width[2]):
                for thick in range(s_thick[0],s_thick[1],s_thick[2]):
                    for vel in range(i_vel[0],i_vel[1],i_vel[2]):
                        for slope in range(i_slope[0],i_slope[1],i_slope[2]):
                            for azimut in range(i_azimut[0],i_azimut[1],i_azimut[2]):
                                for depth in range(i_depth[0],i_depth[1],i_depth[2]):
                                    info = [str(loop)]
                                    """
                                    cmd = grass.make_command("r.impact.tsunami", 
                                                        elevation="dtmpatch",
                                                        lake="lake_depth",
                                                        s_rho=rho,
                                                        s_vol=vol,
                                                        s_width=width, 
                                                        s_thick=thick,
                                                        i_east=1527471,
                                                        i_north=5082144,
                                                        i_vel=vel,
                                                        i_slope=slope,
                                                        i_azimut=azimut,
                                                        output="sens_"+str(loop),
                                                        i_depth=depth,
                                                        flags="agw",
                                                        overwrite=True,
                                                        quiet=False)
                                    print " ".join(cmd)
                                    """
                                    wreport = grass.read_command("r.impact.tsunami", 
                                                        elevation="dtmpatch",
                                                        lake="lake_depth",
                                                        s_rho=rho,
                                                        s_vol=vol,
                                                        s_width=width, 
                                                        s_thick=thick,
                                                        i_east=1527471,
                                                        i_north=5082144,
                                                        i_vel=vel,
                                                        i_slope=slope,
                                                        i_azimut=azimut,
                                                        output=out_prefix+str(loop),
                                                        i_depth=depth,
                                                        flags="agw",
                                                        overwrite=True,
                                                        quiet=False)
                                    wave_stats={}
                                    for l in wreport.split("\n"):
                                        w = l.split("=")
                                        if len(w)==2:
                                            wave_stats[w[0]]=float(w[1])
                                    #print wreport
                                    #print wave_stats.keys()
                                    map_stats = getReport(out_prefix+str(loop))
                                    #reporting
                                    info += [ str(p) for p in [rho,vol,width,thick,vel,slope,azimut,depth] ]

                                    for key in stats:
                                        info.append(str(map_stats[key]))
                                    for key in waves:
                                        info.append(str(wave_stats[key]))
                                    report.append(info)
                                    loop += 1
    #write report
    lines=[]
    lines.append(index[0] + "," + ",".join(pars) + "," + ",".join(stats) + "," + ",".join(waves) + "\n")
    for i in report:
        lines.append(",".join(i) + "\n")
    f = open(out_report, 'w')
    f.writelines(lines)
    f.close()
        
    sys.exit()

if __name__ == "__main__":
    options, flags = grass.parser()
    main()

