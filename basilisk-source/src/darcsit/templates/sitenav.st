<div class="sitenav">
  <fieldset>
    <legend>Site</legend>
    <ul>
      <li><a href="$base$/Front Page">Front page</a></li>
      <li><a href="$base$/_index">All pages</a></li>
      <li><a href="$base$/_activity">Recent activity</a></li>
      $if(feed)$
      <li><a href="$base$/_feed/" type="application/atom+xml" rel="alternate" title="ATOM Feed">Atom feed</a> <img alt="feed icon" src="$base$/img/icons/feed.png"/></li>
      $endif$
      <li><a href="$base$/Help">Help</a></li>
    </ul>
    <form action="/_search" method="get" id="searchform">
      <input type="text" name="patterns" id="patterns"/>
      <input type="submit" value="Search"/>
    </form>
  </fieldset>
  <fieldset>
    <legend>Documentation</legend>
    <ul>
      <li><a href="$base$/Tutorial">Tutorial</a></li>
      <li><a href="$base$/src/INSTALL">Installation</a></li>
      <li><a href="$base$/Basilisk C">Basilisk C</a></li>
      <li><a href="$base$/src/README">Solvers and functions</a></li>
      <li><a href="$base$/src/examples/README">Examples</a></li>
      <li><a href="http://groups.google.com/d/forum/basilisk-fr">User forum</a></li>
      <li><a href="$base$/sandbox/documentation">More documentation</a></li>
    </ul>
  </fieldset>
  <fieldset>
    <legend>Development</legend>
    <ul>
      <li><a href="$base$/src/?history">Recent activity</a></li>
      <li><a href="$base$/src/">Source code</a></li>
      <li><a href="https://hub.darcs.net/basilisk/basilisk/browse/src">Darcs Hub</a></li>
      <li><a href="$base$/src/test/README">Test cases</a></li>
      <li><a href="$base$/sandbox/bugs/README">Bug reports</a></li>
      <li><a href="$base$/src/Contributing">How to contribute</a></li>
      <li><a href="$base$/sandbox/">Play in the sandbox</a></li>
    </ul>
  </fieldset>
</div>
