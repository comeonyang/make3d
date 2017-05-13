#!/usr/bin/perl

$Task = $ARGV[0] or $Task='QuadNew';
$PairListFlag = $ARGV[1] or $PairListFlag = 1;
$Fdir = '/afs/cs/group/reconstruction3d/scratch/TestMultipleImage/';
$File = `ls $Fdir$Task/jpg/*.jpg`;
@File = split(/\.jpg/, $File);
#print @File;
$Path = $Fdir.$Task.'/jpg/';
#print $Path."\n";
if (1)
{
# run all possible pair of images
if ($PairListFlag == 0)
{

	for ($i = 0; $i < @File; ++$i)
	{
		$File[$i] =~ s/$Path//;
		$File[$i] =~ s/\n//;
#		print $File[$i];
	}
	for ($i = 0; $i < @File-1; ++$i)
	{
		for ($j = 0; $j < (@File-1); ++$j)
		{
			if ($i < $j)
			{
				$Img1 = $File[$i];
				$Img2 = $File[$j];
				print $Img1." ";
				print $Img2."\n";
				`qsub -l arch=i686 -v $Img1,$Img2,$Task PairModel.sh`;
			}
		}
	}
}
else
{
# run Pairs im PairList only
	$file = $ARGV[2] or  $file = 'PairList.txt';
	open(INFO, $file);
	@PairList  = <INFO>;
	close(INFO);

	for ( $i = 0; $i < @PairList; ++$i)
	{
		@Target = split(/ /, $PairList[$i]);
		$Target[0] =~ s/\n//;
		$Target[1] =~ s/\n//;
		`export Img1 = $Target[0]`;
		`export Img2 = $Target[1]`;
		`echo '$Img1'`;
		`echo '$Img2'`;
		#print $Img1." ";
		#print $Img2."\n";
		`qsub -l arch=i686 -v '$Img1','$Img2','$Task' PairModel.sh`;
	}
}
}
