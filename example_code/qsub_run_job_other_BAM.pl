@pool_b = (6);
#@pool_l = (50,70);
@pool_l = (70);

$folder = '/dev/ngs005/scratch/lwu/Data/Roche_Batch2/';
$outputs_folder = 'result_Roche_Batch2_BAM/';
opendir(DH,$folder);
@files=readdir(DH);
close DH;
system('mkdir -p '.$outputs_folder);

foreach(@files){
	if (/^(S.*)\.bam$/){
		$file_name = $_;
		$outputs_raw = $1;
		#print $outputs."\n";
		foreach $param_b (@pool_b){
			#print $param_b;
			foreach $param_l (@pool_l){
				$outputs = $outputs_raw.'.b'.$param_b.'.l'.$param_l;
				print $outputs."\n";
				system('qsub SIST_QSUB_BAM.sh '.$folder.
						' '.$file_name.
						' '.$outputs_folder.' '.$outputs.
						' '.$param_b.' '.$param_l);
			}
		}
	}
}
close OFH;