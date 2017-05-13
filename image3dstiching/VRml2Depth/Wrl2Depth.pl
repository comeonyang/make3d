#!/usr/bin/perl

# Read the whole file
$file = $ARGV[0] or  $file = 'building0010';

open(INFO, $file.'.wrl');
@WholeWrl  = <INFO>;
close(INFO);
$IndexedFaceSetStart = 0;
$IndexedFaceSetEnd = 0;
$TextureCoordinateStart = 0;
$TextureCoordinateEnd = 0;

if (1)
{
for ( $i = 0; $i < @WholeWrl; ++$i)
{
	if ( $WholeWrl[$i] =~ /geometry IndexedFaceSet {/)
	{
		print "Find IndexedFaceSetStart"."\n"	;
		$IndexedFaceSetStart = $i+3;
		print @WholeWrl[$IndexedFaceSetStart];
	}
	elsif ( $WholeWrl[$i] =~ /texCoord TextureCoordinate {/)
	{
		print "Find TextureCoordinateStart"."\n"	;
		$TextureCoordinateStart = $i+2;	
		print @WholeWrl[$TextureCoordinateStart];
	}
	elsif ($WholeWrl[$i] =~ /]/)
	{
		if ($IndexedFaceSetEnd == 0 && $IndexedFaceSetStart != 0)
		{
			print "Find End of IndexedFaceSet"."\n";
			$IndexedFaceSetEnd = $i-1;
			print @WholeWrl[$IndexedFaceSetEnd];
		}	
		if ($TextureCoordinateEnd == 0 && $TextureCoordinateStart != 0)
		{
			print "Find End of TextureCoordinate"."\n";
			$TextureCoordinateEnd = $i-1;
			print @WholeWrl[$TextureCoordinateEnd];
		}	
	}
}
@IndexedFaceSet = @WholeWrl[$IndexedFaceSetStart..$IndexedFaceSetEnd];
# clear up the data
for ($i = 0; $i < @IndexedFaceSet; ++$i)
{
	$IndexedFaceSet[$i] =~ s/,//;
}
$OutIndexedFaceSetFile = 'IndexedFaceSet_'.$file;
open(INFO, ">$OutIndexedFaceSetFile");
print INFO @IndexedFaceSet;
close(INFO);

@TextureCoordinate = @WholeWrl[$TextureCoordinateStart..$TextureCoordinateEnd];
for ($i = 0; $i < @TextureCoordinate; ++$i)
{
	$TextureCoordinate[$i] =~ s/,//;
}
$OutTextureCoordinateFile = 'TextureCoordinate_'.$file;
open(INFO, ">$OutTextureCoordinateFile");
print INFO @TextureCoordinate;
close(INFO);

}
