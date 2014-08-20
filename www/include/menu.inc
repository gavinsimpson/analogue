<?php
   echo "
   <div  id=\"menu\">
			<ul>
				<li" . ($current == "home" ? " id=\"current\"" : "") . "><a href=\"index.php\">Home</a></li>
				<li" . ($current == "news" ? " id=\"current\"" : "") . "><a href=\"news.php\">Project News</a></li>
				<li" . ($current == "downloads" ? " id=\"current\"" : "") . "><a href=\"downloads.php\">Downloads</a></li>
				<li" . ($current == "documentation" ? " id=\"current\"" : "") . "><a href=\"documentation.php\">Documentation</a></li>
				<li" . ($current == "development" ? " id=\"current\"" : "") . "><a href=\"devel.php\">Development</a></li>
				<li" . ($current == "support" ? " id=\"current\"" : "") . "><a href=\"support.php\">Support</a></li>	
			</ul>
		</div>
   ";
?>