<?php
header("Cache-Control: no-cache, must-revalidate");
header("Expires: Sat, 26 Jul 1997 05:00:00 GMT");
echo json_encode(array("v" => trim(`md5sum index_raw.html | awk '{print $1}'`)));
?>
