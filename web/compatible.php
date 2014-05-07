<?php
switch($_GET['c']){
	case '0': case '1':
		$c = $_GET['c'];
		echo "Thank you for your help; your browser has been logged as ".($_GET['c']=='1' ? "compatible." : "incompatible - we will work on correcting this. If you would like to provide more information, <a href='https://github.com/aschlosberg/CompressGV/issues'>please lodge an issue</a>.");
		break;
	default;
		exit;
}
$agent = escapeshellarg($_SERVER["HTTP_USER_AGENT"]);
`echo "$c $agent" >> compatible.log`;
?>
