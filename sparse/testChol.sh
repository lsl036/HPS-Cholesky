###
 # @Description: 
 # @Author: Shengle Lin
 # @Date: 2022-05-23 11:16:18
 # @LastEditors: Shengle Lin
 # @LastEditTime: 2022-10-05 18:35:20
### 
IFS_old=$IFS
IFS=$'\n'
for line in $(<alldata.txt)
do
	#echo "$line"
	# arr=($line) 
	# echo "${arr[0]}:${arr[1]}"
	var1=`echo $line|awk -F ' ' '{print $1}'`
	# echo $var1
	var2=`echo $line|awk -F ' ' '{print $2}'`
	# echo $var2
	var3=`echo $line|awk -F ' ' '{print $3}'`
	# var4=`echo $line|awk -F ' ' '{print $4}'`
	# var5=`echo $line|awk -F ' ' '{print $5}'`
	./sparse_choltest $var1 $var2 $var3
done;
IFS=$IFS_old