#!/bin/sh  
# \
exec tclsh "$0" "$@"

puts "this has to be ran from zoot or one of the condor head"

# scan the quad in parallel
# run same file PROC to process the data and generate a 
# file suitable for analysis
set NumI 30
# below S1 is the bunchlength
set S1Min 0.005
set S1Max 0.025

set rootname sigzscan
set index 0
set LoopIndex 0
for { set i 1 } { $i <= $NumI } { incr i } { 
       set Sigma [expr $S1Min + ($S1Max - $S1Min)/($NumI*1.00-1.00)*($i*1.00-1.00)]
       set S1 [expr sqrt($S2/(3.1415*$Sigma))]
       set Sparm $rootname$i.param
       set fs [open $Sparm w]
       puts $fs "$i $Sigma"
       close $fs
       
       puts "SCAN: S1:: $S1 S2:: $S2 ----- $index  Loop $LoopIndex Sigma:: $Sigma"

       set ShFile $rootname$i\_$j.sh
    
       set SCRATCH_Dir /tmp/batchtmp_$rootname$index
       set SOURCE_Dir  [pwd] 
       set fs [open $ShFile w]
       puts $fs "#!/bin/bash"
       puts $fs "# $ShFile"
       puts $fs "source /opt/nicadd/bin/K2setup all"
       puts $fs "export LD_LIBRARY_PATH=\"/usr/lib/openmpi/1.4-gcc/lib/\""
       puts $fs "mkdir $SCRATCH_Dir"
       puts $fs "cd $SCRATCH_Dir"
       puts $fs "cp $SOURCE_Dir/INPUT/* $SCRATCH_Dir/."
       puts $fs "cp $SOURCE_Dir/$rootname$i\_$j * $SCRATCH_Dir/."
       puts $fs "/opt/nicadd/contrib/piot/Argonne_32/epics/extensions/bin/linux-x86/replaceText a0ellipsscan.ref $rootname$i\_$j.in -original=SPOT,QQ -replacement=$S1,$S2"
       puts $fs "/opt/nicadd/contrib/piot/Astra/AstraBetaLSR $rootname$i\_$j.in "


       puts $fs "cp $SCRATCH_Dir/$rootname$i\_$j* $SOURCE_Dir/."
       puts $fs "rm -rf $SCRATCH_Dir"
       close $fs
       exec chmod 775 $ShFile
    
       set tmpFile $rootname$i\_$j.cmd

       set fd [open $tmpFile w]
       puts $fd "# $tmpFile"
       puts $fd "executable                = ./$ShFile"
       puts $fd "requirements              =  Arch==\"X86_64\" "
#       puts $fd "requirements              = (Machine !=\"node1.nicadd.niu.edu\")&&(Machine !=\"node2.nicadd.niu.edu\")&&(Machine !=\"node3.nicadd.niu.edu\")&&(Machine !=\"node4.nicadd.niu.edu\")&&(Machine !=\"node5.nicadd.niu.edu\")&&(Machine !=\"node6.nicadd.niu.edu\")&&(Machine !=\"node7.nicadd.niu.edu\")&&(Machine !=\"node8.nicadd.niu.edu\")&&(Machine !=\"node9.nicadd.niu.edu\")&&(Machine !=\"node10.nicadd.niu.edu\")&&(Machine !=\"node11.nicadd.niu.edu\")&&(Machine !=\"node13.nicadd.niu.edu\")"
       puts $fd "log                       = $rootname$i\_$j.log"
       puts $fd "output                    = $rootname$i\_$j.out"
       puts $fd "error                     = $rootname$i\_$j.error"
       puts $fd "universe                  = vanilla"
       puts $fd "should_transfer_files     = YES"
       puts $fd "when_to_transfer_output   = ON_EXIT"
       puts $fd "notification              = NEVER"
       puts $fd "Queue"
       close $fd

       catch {exec condor_submit $tmpFile } result
       puts "job $index ::  $result"

  }
}
