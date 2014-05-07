<?php
switch($_GET['c']){
	case '0': case '1':
		$c = $_GET['c'];
		echo "Thank you for your help; your browser has been logged as ".($_GET['c']=='1' ? "compatible." : "incompatible - we will work on correcting this. If you would like to provide more information, <a href='https://github.com/aschlosberg/CompressGV/issues/new'>please lodge an issue</a>.");
		break;
	default;
		header("HTTP/1.1 500 Internal Server Error");
		exit;
}
$agent = escapeshellarg($_SERVER["HTTP_USER_AGENT"]);
`echo "$c $agent" >> compatible.log`;
?>
