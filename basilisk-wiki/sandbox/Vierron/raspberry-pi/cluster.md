/** # How to make your own cluster with Raspberry pi ?
![Raspberry-pi cluster](https://drive.google.com/uc?export=view&id=1B9Ds91n-lVa6HuU_ioL7Dmo0bI5DaCev){width=50%}

![Raspberry-pi cluster cases](https://drive.google.com/uc?export=view&id=1J3oTRhsvRHLm_D1sq8HVFv0MYehtkXfE){width=30%}
*/
/** ## What's Raspberry pi ?
[Official website](https://www.raspberrypi.org/help/what-%20is-a-raspberry-pi/)
*/

/** ## Which Raspberry pi ?
Raspberry pi foundation have made a lot of computers during the years [(All products)](https://www.raspberrypi.org/products/).
For our project we will need [Raspberry pi 4B 8Gb](https://www.raspberrypi.org/products/raspberry-pi-4-model-b/)
*/

/** ### Caracteristic of Raspberry pi 4B
   * Broadcom BCM2711, Quad core Cortex-A72 (ARM v8) 64-bit SoC @ 1.5GHz
   * 2GB, 4GB or 8GB LPDDR4-3200 SDRAM (depending on model)
   * 2.4 GHz and 5.0 GHz IEEE 802.11ac wireless, Bluetooth 5.0, BLE
   * Gigabit Ethernet
   * 2 USB 3.0 ports; 2 USB 2.0 ports.
   * Raspberry Pi standard 40 pin GPIO header (fully backwards compatible with previous boards)
   * 2 Ã— micro-HDMI ports (up to 4kp60 supported)
   * 2-lane MIPI DSI display port
   * 2-lane MIPI CSI camera port
   * 4-pole stereo audio and composite video port
   * H.265 (4kp60 decode), H264 (1080p60 decode, 1080p30 encode)
    OpenGL ES 3.0 graphics
   * Micro-SD card slot for loading operating system and data storage
   * 5V DC via USB-C connector (minimum 3A*)
   * 5V DC via GPIO header (minimum 3A*)
   * Power over Ethernet (PoE) enabled (requires separate PoE HAT)
   * Operating temperature: 0 â€“ 50 degrees C ambient
   
   If you are interessted about Raspberry pi you will find a lot of differents projects. [For example](https://pimylifeup.com/category/projects/)
*/

/** ## Hardware
  * At least 2x Raspberry pi
  * Usb Hub charger 10W/port
  * x usb-c charger cables
  * x Micro_SD cards > 16 GB
  * USB drive
  * Optional :
  * * Network switch
  * * x RJ45 cables
  * * Cluster case
  
It's possible to use power over (poe) ethernet on raspberry pi but we need :
Poe-hat + x Poe port Network switch
*/
/** ## Software installation
You can find here a good [installation guide]( https://glmdev.medium.com/building-a-raspberry-pi-cluster-784f0df9afbd)

few commands:

  * get ip address of your Raspberry : 
  
~~~Bash  
sudo nmap -sP 192.168.1.* | grep -A1 'Rasp' | grep '192' | sed -e 's/.*(\(.*\))/\1/'
~~~
*/
/** ## How to launch simulation
[compile your code](http://www.basilisk.fr/src/Tips#running-on-supercomputers)

Now we will need sbatch command to launch our code on raspberry pi 's cluster.
*/
/** ## Overclock cpu
[official guide](https://magpi.raspberrypi.org/articles/how-to-overclock-raspberry-pi-4)

First install small software for stress your cpu

~~~Bash
sudo apt-get install stress
~~~

Small script for stress your cpu and print frequency and temperature of your cpu.

~~~Bash
#!/bin/bash
while (($SECONDS<=1000))
 do
	core1=`cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq`
	core2=`cat /sys/devices/system/cpu/cpu1/cpufreq/scaling_cur_freq`
	core3=`cat /sys/devices/system/cpu/cpu2/cpufreq/scaling_cur_freq`
	core4=`cat /sys/devices/system/cpu/cpu3/cpufreq/scaling_cur_freq`
	temperature=`cat /sys/class/thermal/thermal_zone0/temp`
	echo $SECONDS $core1 $core2 $core3 $core4 $temperature | tee -a data
	sleep 10
done& stress -c 4 -t 900s
~~~
![Graph evolution cpu](https://drive.google.com/uc?export=view&id=119i1V4EqSt2WNBhfyuuwasMB2ubh6RIQ){width=50%}

~~~gnuplot Evolution Frequency and Temperature cpu
#set output 'test-cpu.png'
set grid
set xlabel 'Time (s)'
set ylabel 'Frequency (Mhz)'
set y2tics
set y2label 'Temperature (Â°C)'
plot 'data' u 1:($3/1000000) w l lw 1.5 t 'cpu',\
	'' u 1:($6/1000) w l lw 1.5 axes x1y2 t 'temperature'
~~~

How far can we overclock our cpu ?

It will depend on your cooling solution and your power supply.
The cpu's temperature must be < 80 Â°C and more your cpu's temperature is lower better yout cpu's life time will be.
*/
/** ## Benchmark

I choose [karman.c](http://www.basilisk.fr/src/examples/karman.c) for making my speed test.

The script I used on my computer :

~~~Bash
#!/bin/bash
process_id=$!
CC='mpicc -D_MPI=2' make karman.tst
wait $process_id
echo $SECONDS
~~~

The job script / submission on raspberry pi's cluster :

~~~Bash
#!/bin/bash
#SBATCH -n 2

process_id=$!
mpirun -np 2 ./kar
wait $process_id
echo $SECONDS
~~~

(cat slurm*.out to see the result of echo $SECONDS)

Table: cpu speed inventory 

|CPU |Time (s)  |
|--- | --- | ---|
|i5_6300 2.4GHz| 189 s|
|Broadcom BCM2711, Quad core Cortex-A72 (ARM v8) 64-bit SoC @ **1.5 GHz**|991 s|
|Broadcom BCM2711, Quad core Cortex-A72 (ARM v8) 64-bit SoC @ **2.2 GHz**|797 s|

With ~46.7% increase cpu's speed we save ~20% time calculation.
*/
/** ## "Water cooling" solution
Next step is to put raspberry-pi's cluster in a mineral oil tank (dielectric liquid).
Calorific capacity of mineral oil is twice better than the air.
*/
/** ## Related posts
[Basilisk on Raspberry Pi 4](http://www.basilisk.fr/sandbox/M1EMN/PI4/README)
*/
/** ## To do list
Install 64-bits operating system instead of 32-bits to save some ram.
Like ubuntu arm 64-bits.
*/