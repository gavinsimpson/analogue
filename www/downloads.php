<?php
    // constants to customise this pages headers
    $DESCRIPTION = "Details on downloading the analogue package for R";
    $TITLE = "analogue &mdash; download and install the package";
    $AUTHOR = "Gavin L. Simpson";
    $KEYWORDS = "analogue, R, modern analogue technique, analogue matching, transfer functions, palaeoecology, palaeolimnology";
?>

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<!-- head tags here -->
<?php
    include_once("./include/head.inc");
?>

<body>
<!-- wrap starts here -->
<div id="wrap">
		
		<!--header -->
		<div id="header">			
				
			<h1 id="logo-text">ana<span class="gray">logue</span></h1>		
			<h2 id="slogan">a palaeoecological data analysis package for R</h2>
		</div>
		
		<!-- menu -->		
		<?php
            include_once("./include/menu.inc");
		?>			
			
		<!-- content-wrap starts here -->
		<div id="content-wrap">
				
		<!-- side bar -->
        <?php
            include_once("./include/side_bar.inc");
        ?>
				
			<div id="main">
				
				<a name="Welcome"></a>
				<h1>Downloads</h1>
				
				<h2>Stable Version</h2>
				<p>The current stable release of <strong>analogue</strong> is version 
				<?php echo '{$stable_version}'; ?>. This version is 
				<a href="http://cran.r-project.org/web/packages/analogue/index.html">available from</a> the 
				CRAN mirrors.</p>
				
				<p>To install the stable version of <strong>analogue</strong>, within R type:<br />
				<code>install.packages("analogue", depend = TRUE)</code></p>
				
				<h2>Development Snapshot</h2>
				<p>The development version of analogue is available from the R-Forge repository. The source
				code can be downloaded via SVN. If the development version passes the nightly package checks
				run by the R-Forge system, then packaged source code and binaries for Linux, Windows and
				MacOS X can be installed within R by typing:<br />
				<code>install.packages("analogue", repos="http://R-Forge.R-project.org")</code><br />Note that
				you also need the <strong>vegan</strong> package installed.</p>

			</div>
		
		<!-- content-wrap ends here -->	
		</div>
					
		<!--footer starts here-->
		<?php
            include_once("./include/footer.inc");
        ?>

<!-- wrap ends here -->
</div>

</body>
</html>
