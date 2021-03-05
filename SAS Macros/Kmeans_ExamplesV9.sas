/*EXECUTE THE TRAJECTORIES MACRO FILE BEFORE RUNNING THESE EXAMPLES*/

/*  
Debugging Options -- Generates a very long log!
options mprint symbolgen spool mlogic;run;
*/

/*  Clear Graphics catalog
proc greplay igout = work.gseg nofs;
delete _all_;
run;
quit;
*/
/**************************************************************************************************/


/*Sample Code for Executing a One-Step K-Means Analysis*/
libname traj  "C:\Users\gagnon\Dropbox\Theory\Trajectories\Data generation";
options ls=80 ps=55;

/* Set up graphics destination and type
%macro graphset(gpath,dpi,style,format=emf,type=listing,);*/

%graphset(gpath=C:\Users\gagnon\Documents\figure output,dpi=300,style=splinecurve,
          format=jpeg,type=listing);

/* Initialize trajectories macros
%macro trajsetup(dsn,id,time,riskvar,ngroups,maxdf,ptrim,steps=1,method=0,random=YES);*/

%trajsetup(dsn=traj.Gdata,id=id,time=time,riskvar=response,ngroups=4,maxdf=30,ptrim=0,
           steps=1,method=0,random=YES);
  run;

  /* Do the iterative K-means algorithm
  %macro trajloop(outlib,iter,minchange,minsubs,min_x,max_x,by_x,min_y,max_y,by_y,
                    showall=NO,showany=YES);*/

%trajloop(outlib=traj,iter=50,minchange=0.5,minsubs=0,min_x=-90,max_x=3600,by_x=360,
min_y=-76,max_y=120,by_y=20,showall=NO,showany=YES);
run;
data example1;
do i = 1 to 3;
 call sound(659,200);
 call sound(587,200);
 call sound(523,200);
 end;
run;
proc freq data=

/*Optional macro for generating silhouette plots*/
%silhouette(lib=traj,sid=id,groupsn=4);
run;

/* Macro for finer control of trajectory plots
%macro trajplot(dsn,rvar,groups,timevar,attrid,xlabel,ylabel,min_x,max_x,by_x,min_y,max_y,by_y);*/
%trajplot(graphout,pred,group,age,scatter,Age in Years,SBP,30,100,10,100,200,10);


/*Restarting the system without a random assignment of subjects to clusters, 
use clusout data set created in last iteration and "random=NO" option*/

%graphset(C:\Users\gagnon\Documents\figure output,,splinecurve,format=jpeg,type=listing);

%trajsetup(traj.clusout,newid,age,sbp,6,10,0,steps=1,method=0,random=NO);
  run;
%silhouette(lib=traj,sid=newid,groupsn=6);
run;




/**************************************************************************************************/

/*Sample Code for Executing a Two-Step K-Means/Hierarchical Analysis*/


%graphset(C:\Users\gagnon\Documents\figure output,,splinecurve,format=jpeg,type=listing);
options nomprint nosymbolgen nospool nomlogic;


/*%macro trajsetup(dsn,id,time,riskvar,ngroups,maxdf,ptrim,steps=1);*/
/*options spool mprint mlogic symbolgen;*/
%trajsetup(traj.testing,newid,age,sbp,15,10,0,steps=2,method=0);
  run;
/*%macro trajloop(outlib,iter,minchange,minsubs,min_x,max_x,by_x,min_y,max_y,by_y,
                  showall=NO,showany=YES);*/

%trajloop(traj,50,0.5,10,20,110,20,100,180,20,showany=NO,showall=NO);
run;

/*Hierarchical clustering using trajloop two-step output:*/
/*%macro trajhclus(lib,numclus,res,rfbase);*/
options spool mprint mlogic symbolgen;
options nomprint nosymbolgen nospool nomlogic;

%trajhclus(traj,6,0.1,sbp);


/*%macro trajplot(dsn,rvar,groups,timevar,attrid,xlabel,ylabel,min_x,max_x,by_x,min_y,max_y,by_y);*/

%trajplot(traj.allplot,pred,cluster,age,scatter,Age in Years,SBP,30,100,10,100,200,10);

