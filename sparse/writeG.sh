###
 # @Description: 
 # @Author: Shengle Lin
 # @Date: 2022-05-23 11:16:18
 # @LastEditors: Shengle Lin
 # @LastEditTime: 2022-09-14 12:02:42
### 
IFS_old=$IFS
IFS=$'\n'
for line in $(<GCNdata_write.txt)
do
	#echo "$line"
	# arr=($line) 
	# echo "${arr[0]}:${arr[1]}"
	var1=`echo $line|awk -F ' ' '{print $1}'`
	# echo $var1
	var2=`echo $line|awk -F ' ' '{print $2}'`
	# echo $var2
	var3=`echo $line|awk -F ' ' '{print $3}'`
	var4=`echo $line|awk -F ' ' '{print $4}'`
	var5=`echo $line|awk -F ' ' '{print $5}'`
	./sparse_cholwrite $var1 $var2 $var3 $var4 $var5
done;
IFS=$IFS_old