<form action="$pageUrl$" method="post" id="BogoForm">
  <fieldset>
    <label>Username </label>
    <input type="text" name="bogoid" id="bogoid" size="15" tabindex="1"/>
    <label>Password </label>
    <input type="password" name="bogopid" id="bogopid" size="15"
           tabindex="2"/>
    <input type="text" name="destination" id="destination"
           style="display: none;" value="$referer$"/>
    <input type="submit" name="bogo" id="bogo" value="Login" tabindex="3"/>
  </fieldset>
</form>
<p>If you do not have an account,
  <a href="/_register?destination=$referer$">click here to get one.</a>
</p>
<!-- <p>If you forgot your password, -->
<!--   <a href="/_resetPassword">click here to get a new one.</a></p> -->
