<?php
$files = array("msa", "del", "neut", "classify");
$pref = "CompressGV-".time()."-".substr(md5(uniqid()), 0, 8);
foreach($files as $k => $f){
	$files[$k] = "/tmp/{$pref}_{$f}";
	file_put_contents($files[$k], trim($_POST[$f]));
}

$exec = "./grantham ".implode(" ", $files);
echo json_encode(array("response" => `$exec 2>&1`));
?>
