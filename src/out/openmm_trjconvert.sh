
  # 你的输入文件的路径+文件名
  filename1="./all0.pdb"
  
  
  outname1=md_0
  timename1=time_0.txt
  
  
  echo $filename1  $outname1 $timename1
  echo 0 | gmx trjconv   -f $filename1 -s ../md.tpr -o temp-0.xtc
  echo 0 | gmx trjconv   -f temp-0.xtc -s ../md.tpr -o temp-0.pdb

#  这里面 0  50  0 表示你要设定的pdb文件的开始（0ns）、间隔（50 ps）、开始的frame号（0）
  perl cal_time.pl  temp-0.pdb  temp-1.pdb  0  1  0  
  
  grep "t=" temp-1.pdb > $timename1

  echo 1 0 | gmx trjconv   -f temp-1.pdb  -s ../md.tpr -o temp-2.xtc -pbc cluster 
  echo 1 0 | gmx trjconv   -f temp-2.xtc  -s ../md.tpr -o $outname1.xtc -fit rot+trans

  rm temp*.pdb
  rm temp*.xtc


