@pool_b = (6);
#@pool_l = (50,70);
@pool_l = (70);

$folder = '/storage2/lwu/SEQC2/SIST/Thermo_cfDNA/';
$outputs_folder = 'result_Thermo_cfDNA/';
opendir(DH,$folder);
@files=readdir(DH);
close DH;
system('mkdir -p '.$outputs_folder);

foreach(@files){
	if (/^[A-Z]/){
		my $sub_folder = $_;
		opendir(DH_sub, $folder.'/'.$sub_folder);
		@files_sub = readdir(DH_sub);
		close DH_sub;
		
		foreach(@files_sub){
			if (/^(Tag.*)\.bam$/){
				$file_name = $_;
				$outputs_raw = $1;
				#print $outputs."\n";
				foreach $param_b (@pool_b){
					#print $param_b;
					foreach $param_l (@pool_l){
						$outputs = $sub_folder.'.'.$outputs_raw.'.b'.$param_b.'.l'.$param_l;
						print $outputs."\n";
						system('qsub SIST_QSUB_BAM.sh '.$folder.
								' '.$sub_folder.'/'.$file_name.
								' '.$outputs_folder.' '.$outputs.
								' '.$param_b.' '.$param_l);
					}
				} 
			}
		}
	}
}
close OFH;