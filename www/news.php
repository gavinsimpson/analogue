<?php
    // constants to customise this pages headers
    $DESCRIPTION = "News on the analogue project";
    $TITLE = "analogue &mdash; Project News";
    $AUTHOR = "Gavin L. Simpson";
    $KEYWORDS = "analogue, R, modern analogue technique, analogue matching, transfer functions, palaeoecology, palaeolimnology";
    $current = "news";
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
				<h1>analogue News</h1>
			    <h2>Version 0.6-8 uploaded to CRAN</h2>
			    <h3>7 May 2009</h3>
			    <p>A new version of <strong>analogue</strong> was released in order to fix a problem with the CITATION
			    on CRAN. Several additional new functions and methods were included in this release as a result:</p>
			    <ul>
        	        <li><strong><kbd>residLen</kbd></strong>&mdash;new function to compute the squared residual length diagnostic
        	        for passive samples in a constrained ordination. Includes several plotting routines for both base
        	        Lattice graphics.</li>
        	        <li><strong><kbd>stdError</kbd></strong>&mdash;computes the weighted standard deviation of the predicted values
        	        over the <em>k</em>-closest analogues. This measure has been proposed as an uncertainty measure for
        	        MAT model predictions. Methods provided for <strong><kbd>mat</kbd></strong> and 
        	        <strong><kbd>predict.mat</kbd></strong></li>.
        	        <li><strong><kbd>predict.mat</kbd></strong> now returns the dissimilarity matrix between training 
        	        set and new samples.</li>
        	        <li><strong><kbd>getK</kbd></strong>&mdash;new method for <strong><kbd>predict.mat</kbd></strong>.</li>
			    </ul>
			    <p>See the <a href="http://cran.r-project.org/web/packages/analogue/ChangeLog">ChangeLog</a>
			    for more details.</p>
				<h2>Version 0.6-6 uploaded to CRAN</h2>
				<h3>10 March 2009</h3>
				<p>A new version of analogue has been uploaded to CRAN. It contains lots of new functionality 
				and fixes, especially to the ROC curve functions. See the ChangeLog for details of the fixes and 
				improvements.  For version 0.7-0 I hope to tidy the package and help up a bit, add extra documentation 
				and vignettes, and add a proper NEWS file.</p>

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