<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
          "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    $if(basepage)$
    <base href="$base$/$basepage$/"/>
    $endif$
    <title>$wikititle$ - $pagetitle$</title>
    <link href="$base$/css/custom.css" rel="stylesheet" media="screen, projection" type="text/css" />
    <link href="$base$/css/print.css" rel="stylesheet" media="print" type= "text/css" />
    <!--[if IE]><link href="$base$/css/ie.css" rel="stylesheet" media="screen, projection" type="text/css" /><![endif]-->
    <link rel="stylesheet" href="$base$/css/basilisk.css"/>
    <script src="$base$/js/jquery.min.js" type="text/javascript"></script>
    <script src="$base$/js/jquery-ui.packed.js" type="text/javascript"></script>
    <script src="$base$/js/plots.js" type="text/javascript"></script>
    $math()$
  </head>
  <body>
    <div id="doc3" class="yui-t1">
      <div id="yui-main">
        <div id="maincol">
	  <div id="content">
	    $if(notitle)$
	    $else$
	    <h1 class="pageTitle"><a href="$pageUrl$">$pagetitle$</a></h1>
	    $endif$
	    <link rel="stylesheet" href="$base$/css/edit.css"/>
	    <link rel="stylesheet" href="$base$/css/codemirror.css"/>
	    <link rel="stylesheet" href="$base$/codemirror/addon/display/fullscreen.css"/>
	    <link rel="stylesheet" href="$base$/codemirror/addon/fold/foldgutter.css" />
	    <script src="$base$/codemirror/lib/codemirror.js"></script>
	    <script src="$base$/codemirror/mode/clike/clike.js"></script>
	    <script src="$base$/codemirror/mode/markdown/markdown.js"></script>
	    <script src="$base$/codemirror/addon/display/fullscreen.js"></script>
	    <script src="$base$/codemirror/addon/edit/matchbrackets.js"></script>
	    <script src="$base$/codemirror/keymap/emacs.js"></script>
	    <script src="$base$/codemirror/addon/comment/comment.js"></script>
	    <script src="$base$/codemirror/addon/dialog/dialog.js"></script>
	    <script src="$base$/codemirror/addon/search/searchcursor.js"></script>
	    <script src="$base$/codemirror/addon/search/search.js"></script>
	    <script src="$base$/codemirror/addon/fold/foldcode.js"></script>
	    <script src="$base$/codemirror/addon/fold/foldgutter.js"></script>
	    <script src="$base$/codemirror/addon/fold/brace-fold.js"></script>
	    <script src="$base$/codemirror/addon/fold/xml-fold.js"></script>
	    <script src="$base$/codemirror/addon/fold/markdown-fold.js"></script>
	    <script src="$base$/codemirror/addon/fold/comment-fold.js"></script>
	    <form action="$pageUrl$" method="post" id="editform"
		  enctype="application/x-www-form-urlencoded">
	      <div id="editcontainer">
		<div id="left-column">
		  <div id="preview">
		    <div id="previewpane"></div>
		  </div>
		  <div id="status"></div>
		</div>
		<div id="right-column">
		  <div id="editor">
		    <textarea cols="80" name="editedText" id="editedText">
		      $if(edit)$
		      $edit$
		      $endif$
		    </textarea>
		  </div>
		  <div id="changes-box">
		    $messages()$
		    <label>Description of changes:</label><br/>
		    <input type="text" name="logMsg" id="logMsg" value="$logMsg$"/><br/>
		    $checkboxes$
		    <input type="button" class="editButton" id="update"
  			   onClick="save();" value="Save"/>
		    &nbsp;<input type="submit" name="cancel" id="cancel" value="Discard"/>
		    &nbsp;<input type="button" class="editButton" id="previewButton"
				 onClick="updatePreviewPane();" style="display: none;"
				 value="Preview"/>
		    &nbsp;<input type="button" class="editButton" id="runButton"
				 onClick="run();" style="display: none;"
				 value="Run"/>
		  </div>
		</div>
	      </div>
	    </form>
	    <script type="text/javascript">
	      var editor;
	      function setup_codemirror() {
	      CodeMirror.keyMap.emacs.Enter = "newlineAndIndent";
	      CodeMirror.extendMode("clike", {electricChars: "{};"});
	      editor = 
    	      CodeMirror.fromTextArea(document.getElementById("editedText"), {
	        matchBrackets: true,
	        mode: "text/x-basiliskc",
		theme: "basilisk",
		tabSize: 8,
		electricChars: true,
		lineWrapping: true,
		foldGutter: true,
		keyMap: "emacs",
      		extraKeys: {
		    "F11": function(cm) {
			cm.setOption("fullScreen", !cm.getOption("fullScreen"));
		    },
		    "Esc": function(cm) {
			if (cm.getOption("fullScreen")) cm.setOption("fullScreen", false);
		    },
		    "Ctrl-E": function(cm){ cm.foldCode(cm.getCursor()); }
		},
		gutters: ["CodeMirror-foldgutter"]
              });
	      editor.on("blur", function() {editor.save()});
	      update();
	      var pending;
	      editor.on("change", function() {
	      clearTimeout(pending);
	      pending = setTimeout(update, 400);
	      });
	      function update() {
	      editor.setOption("mode", /.*\.[ch]/.test("$pageUrl$") ? 
	                       "text/x-basiliskc" : "markdown");
	      }
	      }
	    </script>
	  </div>
	  $body$
	</div>
	$footer()$	  
      </div>
    </div>
    $javascripts$
  </body>
</html>
