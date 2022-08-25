#!/bin/bash

echo

realpath() {
  # $1 : relative filename
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

case=$1
raddir="$(realpath $2)"
tempdir="$(realpath $3)"
outparentdir="$(realpath $4)"
allfigs=$5
simdir=$6
doppdir=$7

ipoldir="$(realpath )"

if [[ ! $case == 2* ]]; then

    configdir=$ipoldir/configtxt/$case

    case $case in
        olym1)
            stt='20151031_0300'
            edt='20151101_0600' ;;
        olym2)
            stt='20151112_1500'
            edt='20151114_0000' ;;
        olym3)
            stt='20151116_1500'
            edt='20151118_0000' ;;
        olym4)
            stt='20151203_0300'
            edt='20151204_0300' ;;
        olym5)
            stt='20151208_0300'
            edt='20151209_1200' ;;
        swbc1)
            stt='20211113_2100'
            edt='20211115_2100' ;;
    esac

else

    st=$(echo $case | cut -d ',' -f1 | sed 's/[^0-9]*//g')
	stdate=${st:0:8}
	sttime=${st:8:2}
    en=$(echo $case | cut -d ',' -f2 | sed 's/[^0-9]*//g')
	endate=${en:0:8}
	entime=${en:8:2}

	stt=${stdate}_${sttime}00
	edt=${endate}_${entime}00

    configdir=$ipoldir/configtxt/${stt}_${edt}

fi

mkdir -p $configdir

echo Selecting radar files for analysis in range $stt to $edt...
echo
sleep 3

station=$(basename $raddir)
tfile=input_${stt}_${edt}.txt

if [[ "$station" == "CASAG" ]]; then
    agency='cwr'
elif [[ "$station" == "NPOL" ]]; then
    agency='olympex'
else
    agency='nexrad'
fi

for filepath in $(ls $raddir/*); do
    file=$(basename $filepath)
    filedt=$(echo $file | cut -d '_' -f2)$(echo $file | cut -d '_' -f3)
    if [ "$filedt" -ge "$(echo $stt | tr -d '_')00" ] && [ "$filedt" -le "$(echo $edt | tr -d '_')00" ]; then 
        echo $filepath >> $configdir/$tfile
        echo $(basename $filepath)
    fi
    if [ "$filedt" -gt "$(echo $edt | tr -d '_')00" ]; then
        break
    fi
done

if [[ "$(head -n 1 $configdir/$tfile | xargs basename)" == "wrfout"* ]]; then
	data='wrf'
else
	data='obs'
fi
inputfile=input_${data}_${stt}_${edt}.txt
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

tfile=temp_${stt}_${edt}.txt

for filepath in $(ls $tempdir/*); do
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
	temp='wrf'
    tempfile=temp_${temp}_mp$(basename $tempdir)_${stt}_${edt}.txt
    snd_on='False'
    wrft_on='True'
else
	temp='snd'
    tempfile=temp_${temp}_${stt}_${edt}.txt
    snd_on='True'
    wrft_on='False'
fi
mv $configdir/$tfile $configdir/$tempfile

if [ -z $simdir ]; then
    fold=$data
else
    fold='obswrf'
fi

if [ -z $doppdir ]; then
    dd_on='False'
else
    dd_on='True'
fi

outdir=$outparentdir/${fold}_temp${temp}_${station}_${stt}_${edt}
mkdir -p $outdir

echo
echo Creating config file for iPOLARRIS...
echo
sleep 3

template=$ipoldir/${agency}_${data}_config_template.txt
configfile=config_${data}_${stt}_${edt}.txt
cp $template $configdir/$configfile

sed -i '' "s/^type ==.*/type == $data == # Type of input data: 'obs' OR 'wrf' (obs + simulated)/g" $configdir/$configfile
sed -i '' "s/.*mphys ==.*/mphys == $data == # Type of microphysics used in model: 'obs' OR '<scheme>' if type = 'wrf'/g" $configdir/$configfile
sed -i '' "s/.*sdatetime ==.*/sdatetime == '$(echo $stt | tr '_' '-')' == # Start time of analysis of interest/g" $configdir/$configfile
sed -i '' "s/.*edatetime ==.*/edatetime == '$(echo $edt | tr '_' '-')' == # End time of analysis of interest/g" $configdir/$configfile
sed -i '' "s%.*rfiles ==.*%rfiles == '$configdir/$inputfile' == # Path to list of radar files to read in%g" $configdir/$configfile
sed -i '' "s%.*wfiles ==.*%wfiles == '$configdir/$tempfile' == # Path to list of WRF temperature files to read in%g" $configdir/$configfile
sed -i '' "s/.*exper ==.*/exper == $station == # Radar location/g" $configdir/$configfile
sed -i '' "s/.*lat ==.*/lat == $latcen == # Latitude of the radar station/g" $configdir/$configfile
sed -i '' "s/.*lon ==.*/lon == $loncen == # Longitude of the radar station/g" $configdir/$configfile
sed -i '' "s%.*image_dir ==.*%image_dir == '$outdir/' == # Output figure directory%g" $configdir/$configfile
sed -i '' "s/.*dd_on ==.*/dd_on == $dd_on == # Doppler gridded velocity on/g" $configdir/$configfile
sed -i '' "s/.*snd_on ==.*/snd_on == $snd_on == # Sounding temperature on/g" $configdir/$configfile
sed -i '' "s/.*wrft_on ==.*/wrft_on == $wrft_on == # WRF temperature on/g" $configdir/$configfile

echo Running iPOLARRIS...
sleep 3

python run_ipolarris.py $configdir/$configfile
