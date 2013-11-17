<?php
/*
Copyright 2013 Arran Schlosberg.

This file is part of https://github.com/aschlosberg/CompressGV (CompressGV)

    CompressGV is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CompressGV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CompressGV. If not, see <http://www.gnu.org/licenses/>.
*/

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
$output = `$exec 2>&1`;

if(!isset($_GET['format'])){
	$_GET['format'] = "json";
}
switch(strtolower($_GET['format'])){
	case "text":
		echo $output;
		break;
	case "json": default:
		echo json_encode(array("response" => $output));
		break;
}

foreach($files as $f){
	unlink($f);
}
?>
