<?php
    echo "
<head>

<meta name=\"Description\" content=\"{$DESCRIPTION}\" />
<meta name=\"Keywords\" content=\"{$KEYWORDS}\" />
<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />
<meta name=\"Distribution\" content=\"Global\" />
<meta name=\"author\" content=\"{$AUTHOR}\"/>
<meta name=\"Robots\" content=\"index,follow\" />

<link rel=\"stylesheet\" href=\"images/Refresh.css\" type=\"text/css\" />

<title>{$TITLE}</title>

</head>
<script type=\"text/javascript\">
var gaJsHost = ((\"https:\" == document.location.protocol) ? \"https://ssl.\" : \"http://www.\");
document.write(unescape(\"%3Cscript src='\" + gaJsHost + \"google-analytics.com/ga.js' type='text/javascript'%3E%3C/script%3E\"));
</script>
<script type=\"text/javascript\">
try {
var pageTracker = _gat._getTracker(\"UA-7900310-1\");
pageTracker._trackPageview();
} catch(err) {}</script>
    ";
?>