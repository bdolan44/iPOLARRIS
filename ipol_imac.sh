#!/bin/bash

echo

realpath() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

ipoldir="$(realpath )"

ptype=$1
starg=$2
enarg=$3
raddir="$(realpath $4)"
tempdir="$(realpath $5)"
simdir=$6
doppdir=$7

mkdir -p $raddir $tempdir
if [ ! -z $simdir ]; then
    mkdir -p $simdir
fi
if [ ! -z $doppdir ]; then
    mkdir -p $doppdir
fi

declare -a mpopts=( "mp06" "mp08" "mp10" "mp16" "mp51" )
declare -a mpnames=( "wsm6" "thom" "morr" "wdm6" "p3" )

st=$(echo $starg | sed 's/[^0-9]*//g')
stdate=${st:0:8}
sttime=${st:8:4}
en=$(echo $enarg | sed 's/[^0-9]*//g')
endate=${en:0:8}
entime=${en:8:4}
stt=${stdate}_${sttime}
edt=${endate}_${entime}

configdir=$ipoldir/configtxt/${stt}_${edt}
mkdir -p $configdir

echo Determining data agency and type...
echo
sleep 3

station=$(basename $raddir)
if [[ "$station" == "CASAG" ]]; then
    agency='cwr'
elif [[ "$station" == "NPOL" ]]; then
    agency='olympex'
else
    agency='nexrad'
fi
echo $agency

if [[ "$(ls $raddir/* | head -n 1 | xargs basename)" == "wrfout"* ]]; then
	data='wrf'
    mp=$(ls $raddir/* | head -n 1 | xargs basename | cut -d '_' -f2)
    for ((ii=0;ii<${#mpopts[@]};ii++)); do
       [[ "${mpopts[ii]}" = "$mp" ]] && break
    done
    inputfile=input_${mpnames[ii]}_${stt}_${edt}.txt
    configfile=config_${mpnames[ii]}_${stt}_${edt}.txt
    echo $data
else
    if [ -z $simdir ]; then
        data='obs'
        inputfile=input_${data}_${stt}_${edt}.txt
        configfile=config_${data}_${stt}_${edt}.txt
        echo $data
    else
        inputfile=input_obs_${stt}_${edt}.txt
        configfile=config_obs_${stt}_${edt}.txt
        echo "obs vs wrf"
    fi
fi

echo
echo Selecting radar files for analysis in range $stt to $edt...
echo
sleep 3

tfile=input_${stt}_${edt}.txt

for filepath in $(ls $raddir/* | sort); do
    file=$(basename $filepath)
    if [[ "$(ls $raddir/* | head -n 1 | xargs basename)" == "wrfout"* ]]; then
        filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
    else
        filedt=$(echo $file | cut -d '_' -f2)$(echo $file | cut -d '_' -f3)
    fi
    if [ "$filedt" -ge "$(echo $stt | tr -d '_')00" ] && [ "$filedt" -lt "$(echo $edt | tr -d '_')00" ]; then 
        echo $filepath >> $configdir/$tfile
        echo $(basename $filepath)
    fi
    if [ "$filedt" -ge "$(echo $edt | tr -d '_')00" ]; then
        break
    fi
done

mv $configdir/$tfile $configdir/$inputfile

#latcen=$(ncdump $indir/$headfile | grep "latitude =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)
#loncen=$(ncdump $indir/$headfile | grep "longitude =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)
latcen=47.116943359375
loncen=-124.106666564941
#latcen=49.0164413452148
#loncen=-122.487358093262

echo
echo Selecting temperature files for analysis in range $stt to $edt...
echo
sleep 3

tfile=tmp_${stt}_${edt}.txt

for filepath in $(ls $tempdir/* | sort); do
    file=$(basename $filepath)
    filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
    if [ "$filedt" -ge "$(echo $stt | tr -d '_')00" ] && [ "$filedt" -le "$(echo $edt | tr -d '_')00" ]; then 
        echo $filepath >> $configdir/$tfile
        echo $(basename $filepath)
    fi
    if [ "$filedt" -gt "$(echo $edt | tr -d '_')00" ]; then
        break
    fi
done

if [[ "$(head -n 1 $configdir/$tfile | xargs basename)" == "wrfout"* ]]; then
    mp=$(head -n 1 $configdir/$tfile | xargs basename | cut -d '_' -f2)
    for ((ii=0;ii<${#mpopts[@]};ii++)); do
       [[ "${mpopts[ii]}" = "$mp" ]] && break
    done
    mpname=$(echo "${mpnames[ii]}")
    tempsrc=$mpname
    tempfile=temp_${tempsrc}_${stt}_${edt}.txt
    snd_on='False'
    wrft_on='True'
else
    tempsrc='uwyo'
    tempfile=temp_${tempsrc}_${stt}_${edt}.txt
    snd_on='True'
    wrft_on='False'
fi
mv $configdir/$tfile $configdir/$tempfile

if [ -z $simdir ]; then
    if [[ "$data" == "obs" ]]; then
        fold="obs"
    else
        fold=$mpname
    fi

else
    
    echo
    echo Selecting wrfout files for analysis in range $stt to $edt...
    echo
    sleep 3

    fold='obsvwrf'
    mp2=$(ls $simdir/* | head -n 1 | xargs basename | cut -d '_' -f2)
    for ((ii=0;ii<${#mpopts[@]};ii++)); do
       [[ "${mpopts[ii]}" = "$mp2" ]] && break
    done
    mpname2=$(echo ${mpnames[ii]})

    inputfile2=input_${mpname2}_${stt}_${edt}.txt
    configfile2=config_${mpname2}_${stt}_${edt}.txt
    
    tfile=input_${stt}_${edt}.txt

    for filepath in $(ls $simdir/* | sort); do
        file=$(basename $filepath)
        filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
        if [ "$filedt" -ge "$(echo $stt | tr -d '_')00" ] && [ "$filedt" -lt "$(echo $edt | tr -d '_')00" ]; then 
            echo $filepath >> $configdir/$tfile
            echo $(basename $filepath)
        fi
        if [ "$filedt" -ge "$(echo $edt | tr -d '_')00" ]; then
            break
        fi
    done

    mv $configdir/$tfile $configdir/$inputfile2

fi

if [ -z $doppdir ]; then
    dd_on='False'
else
    dd_on='True'
fi

outfigdir=outputfig/${fold}_temp${tempsrc}_${station}_${stt}_${edt}
outrrdir=$(cd $raddir/../../ && pwd)/radar_rainrates/$station
mkdir -p $outfigdir $outrrdir

if [ -z $simdir ]; then

    echo
    echo Creating config file for iPOLARRIS...
    echo
    sleep 3

    template=$ipoldir/${agency}_${data}_config.txt
    cp $template $configdir/$configfile

    sed -i '' "s/^type ==.*/type == $data == # Type of input data: 'obs' OR 'wrf' (obs + simulated)/g" $configdir/$configfile
    if [[ "$data" == "obs" ]]; then
        sed -i '' "s/.*mphys ==.*/mphys == $data == # Type of microphysics used in model: 'obs' OR '<scheme>' if type = 'wrf'/g" $configdir/$configfile
    else
        sed -i '' "s/.*mphys ==.*/mphys == $mpname == # Type of microphysics used in model: 'obs' OR '<scheme>' if type = 'wrf'/g" $configdir/$configfile
    fi
    sed -i '' "s/.*ptype ==.*/ptype == '$ptype' == # Output figure file extenstion (i.e. png, jpg, mp4, ...)/g" $configdir/$configfile
    sed -i '' "s/.*sdatetime ==.*/sdatetime == '$(echo $stt | tr '_' '-')' == # Start time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s/.*edatetime ==.*/edatetime == '$(echo $edt | tr '_' '-')' == # End time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s%.*rfiles ==.*%rfiles == '$configdir/$inputfile' == # Path to list of radar files to read in%g" $configdir/$configfile
    sed -i '' "s%.*wfiles ==.*%wfiles == '$configdir/$tempfile' == # Path to list of WRF temperature files to read in%g" $configdir/$configfile
    if [[ "$data" == "obs" ]]; then
        sed -i '' "s/.*exper ==.*/exper == $station == # Radar location/g" $configdir/$configfile
    else
        sed -i '' "s/.*exper ==.*/exper == $station-$(echo $mpname | tr '[:lower:]' '[:upper:]') == # Radar location/g" $configdir/$configfile
    fi
    sed -i '' "s/lat ==  == #.*/lat == $latcen == # Latitude of the radar station/g" $configdir/$configfile
    sed -i '' "s/lon ==  == #.*/lon == $loncen == # Longitude of the radar station/g" $configdir/$configfile
    sed -i '' "s%.*image_dir ==.*%image_dir == '$outfigdir/' == # Output figure directory%g" $configdir/$configfile
    sed -i '' "s%.*rr_dir ==.*%rr_dir == '$outrrdir/' == # Output rain rate netcdf directory%g" $configdir/$configfile
    sed -i '' "s/.*dd_on ==.*/dd_on == $dd_on == # Doppler gridded velocity on/g" $configdir/$configfile
    sed -i '' "s/.*snd_on ==.*/snd_on == $snd_on == # Sounding temperature on/g" $configdir/$configfile
    sed -i '' "s/.*wrft_on ==.*/wrft_on == $wrft_on == # WRF temperature on/g" $configdir/$configfile

    echo Running iPOLARRIS...
    sleep 3

    python run_ipolarris.py $configdir/$configfile

else

    echo
    echo Creating OBS and SIM config files for iPOLARRIS...
    echo
    sleep 3

    template=$ipoldir/${agency}_obs_config.txt
    cp $template $configdir/$configfile

    sed -i '' "s/^type ==.*/type == obs == # Type of input data: 'obs' OR 'wrf' (obs + simulated)/g" $configdir/$configfile
    sed -i '' "s/.*mphys ==.*/mphys == obs == # Type of microphysics used in model: 'obs' OR '<scheme>' if type = 'wrf'/g" $configdir/$configfile
    sed -i '' "s/.*ptype ==.*/ptype == '$ptype' == # Output figure file extenstion (i.e. png, jpg, mp4, ...)/g" $configdir/$configfile
    sed -i '' "s/.*sdatetime ==.*/sdatetime == '$(echo $stt | tr '_' '-')' == # Start time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s/.*edatetime ==.*/edatetime == '$(echo $edt | tr '_' '-')' == # End time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s%.*rfiles ==.*%rfiles == '$configdir/$inputfile' == # Path to list of radar files to read in%g" $configdir/$configfile
    sed -i '' "s%.*wfiles ==.*%wfiles == '$configdir/$tempfile' == # Path to list of WRF temperature files to read in%g" $configdir/$configfile
    sed -i '' "s/.*exper ==.*/exper == $station == # Radar location/g" $configdir/$configfile
    sed -i '' "s/lat ==  == #.*/lat == $latcen == # Latitude of the radar station/g" $configdir/$configfile
    sed -i '' "s/lon ==  == #.*/lon == $loncen == # Longitude of the radar station/g" $configdir/$configfile
    sed -i '' "s%.*image_dir ==.*%image_dir == '$outfigdir/' == # Output figure directory%g" $configdir/$configfile
    sed -i '' "s%.*rr_dir ==.*%rr_dir == '$outrrdir/' == # Output rain rate netcdf directory%g" $configdir/$configfile
    sed -i '' "s/.*dd_on ==.*/dd_on == $dd_on == # Doppler gridded velocity on/g" $configdir/$configfile
    sed -i '' "s/.*snd_on ==.*/snd_on == $snd_on == # Sounding temperature on/g" $configdir/$configfile
    sed -i '' "s/.*wrft_on ==.*/wrft_on == $wrft_on == # WRF temperature on/g" $configdir/$configfile

    template2=$ipoldir/${agency}_wrf_config.txt
    cp $template2 $configdir/$configfile2

    sed -i '' "s/^type ==.*/type == wrf == # Type of input data: 'obs' OR 'wrf' (obs + simulated)/g" $configdir/$configfile2
    sed -i '' "s/.*mphys ==.*/mphys == $mpname2 == # Type of microphysics used in model: 'obs' OR '<scheme>' if type = 'wrf'/g" $configdir/$configfile2
    sed -i '' "s/.*ptype ==.*/ptype == '$ptype' == # Output figure file extenstion (i.e. png, jpg, mp4, ...)/g" $configdir/$configfile2
    sed -i '' "s/.*sdatetime ==.*/sdatetime == '$(echo $stt | tr '_' '-')' == # Start time of analysis of interest/g" $configdir/$configfile2
    sed -i '' "s/.*edatetime ==.*/edatetime == '$(echo $edt | tr '_' '-')' == # End time of analysis of interest/g" $configdir/$configfile2
    sed -i '' "s%.*rfiles ==.*%rfiles == '$configdir/$inputfile2' == # Path to list of radar files to read in%g" $configdir/$configfile2
    sed -i '' "s%.*wfiles ==.*%wfiles == '$configdir/$tempfile' == # Path to list of WRF temperature files to read in%g" $configdir/$configfile2
    sed -i '' "s/.*exper ==.*/exper == $station-$(echo $mpname2 | tr '[:lower:]' '[:upper:]') == # Radar location/g" $configdir/$configfile2
    sed -i '' "s/lat ==  == #.*/lat == $latcen == # Latitude of the radar station/g" $configdir/$configfile2
    sed -i '' "s/lon ==  == #.*/lon == $loncen == # Longitude of the radar station/g" $configdir/$configfile2
    sed -i '' "s%.*image_dir ==.*%image_dir == '$outfigdir/' == # Output figure directory%g" $configdir/$configfile2
    sed -i '' "s%.*rr_dir ==.*%rr_dir == '$outrrdir/' == # Output rain rate netcdf directory%g" $configdir/$configfile2
    sed -i '' "s/.*dd_on ==.*/dd_on == $dd_on == # Doppler gridded velocity on/g" $configdir/$configfile2
    sed -i '' "s/.*snd_on ==.*/snd_on == $snd_on == # Sounding temperature on/g" $configdir/$configfile2
    sed -i '' "s/.*wrft_on ==.*/wrft_on == $wrft_on == # WRF temperature on/g" $configdir/$configfile2

    echo Running iPOLARRIS...
    sleep 3

    python run_ipolarris.py $configdir/$configfile $configdir/$configfile2

fi
