<?php
$files = array("msa", "del", "neut", "novel");
$pref = "CompressGV-".time()."-".substr(md5(uniqid()), 0, 8);

$_POST['msa'] = str_replace('&gt;', '>', $_POST['msa']);

foreach($files as $k => $f){
	$files[$k] = "/tmp/{$pref}_{$f}";
	file_put_contents($files[$k], trim($_POST[$f]));
}

if(strpos($_POST['msa'], ">")!==false){
	`./fastaToSingleLine.sh $files[0] > {$files[0]}_single; mv {$files[0]}_single $files[0];`;
}

$exec = "./grantham ".implode(" ", $files);
echo json_encode(array("response" => `$exec 2>&1`));
?>
