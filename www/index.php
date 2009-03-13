<?php
    // constants to customise this pages headers
    $DESCRIPTION = "analogue R package for quantiative palaeoecology, analogue methods and transfer functions";
    $TITLE = "analogue &mdash; an R package for quantitative palaeoecology";
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
				<h1>Welcome to the analogue project</h1>
				
				<p><strong>analogue</strong> is an <a href="http://www.r-project.org">R</a> package for analysing
				palaeoecological data. analogue was started by 
				<a href="http://www.ecrc.ucl.ac.uk/content/view/90/94/">Gavin Simpson</a> to bring together
				functions written for various research projects and his Ph.D thesis. Recently, 
				<a href="http://cc.oulu.fi/~jarioksa/">Jari Oksanen</a> has joined the development team.</p>
				<p><strong>analogue</strong> is free and open-source, released under version 2 of the
				<a href="http://www.gnu.org/licenses/old-licenses/gpl-2.0.html#SEC1">GNU GPL</a>.</p>
				
				<p>The main features of the current version of <strong>analogue</strong> are:</p>
				
				<ul>
            	    <li>Palaeoecological transfer functions:
            	        <ul>
            	            <li>Modern Analogue Technique (<acronym title="Modern Analogue Technique">MAT</acronym>)</li>
            	            <li>Weighted Averaging (<acronym title="Weighted Averaging">WA</acronym>)
            	                <ul>
            	                    <li>Tolerance down-weighting</li>
            	                    <li>Classical and inverse deshrinking, plus others</li>
            	                </ul>
            	            </li>
            	            <li>Transfer function evaluation:
            	                <ul>
            	                    <li>Bootstrap, k-fold, leave-one-out cross validation</li>
            	                    <li>Analogue statistics</li>
            	                </ul>
            	            </li>
            	        </ul>
            	    </li>
            	    <li>Analogue matching</li>
            	    <li>Dissimilarity coefficients, including:
            	        <ul>
                            <li>Chord distance</li>
                            <li>Bray-Curtis distance</li>
                            <li>Gower's General (dis)similarity</li>
                            <li>Manhattan metric</li>
                            <li>...and many more</li>
            	        </ul>
            	    </li>
            	    <li>Methods to select dissimilaity decision thresholds:
            	        <ul>
                            <li>Receiver Operator Characteristic (<acronym title="Receiver Operator Characteristic">ROC</acronym>) curves </li>
                            <li>Monte Carlo resampling</li>
                            <li>Logistic regression modelling</li>            	        
            	        </ul>
            	    </li>
            	    <li>Utilities:
            	        <ul>
            	            <li>Merge data sets</li>
            	            <li>Apply common data transformations</li>
            	        </ul>
            	    </li>
				</ul>
				
				<p>Plenty more features are planned for future versions.</p>

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