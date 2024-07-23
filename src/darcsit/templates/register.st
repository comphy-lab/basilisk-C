Please note that:
<ul>
  <li>Usernames can contain spaces. The prefered format is "Jane Smith".</li>
  <li>Usernames will appear in the file history.</li>
  <li>Email is optional but is recommended for validating authorship.</li>
</ul>
<form action="" method="post" id="bogoForm">
  <fieldset>
    <label>Username (at least 3 letters or digits):</label><br/>
    <input type="text" name="bogoid" id="bogoid"
	   size="20" tabindex="2"/><br/>
    <label>Email (optional, will not be displayed on the Wiki):</label><br/>
    <input type="text" name="email" id="email"
	   size="20" tabindex="3" value=""/><br/>
    <input type="text" name="full_name_1" id="full_name_1"
	   size="20" class="req"/>
    <label>Password (at least 6 characters, including at least one non-letter):</label><br/>
    <input type="password" name="bogopid" id="bogopid"
	   size="20" tabindex="4"/> <br/>
    <label>Confirm Password:</label><br/>
    <input type="password" name="bogopid2" id="bogopid2"
	   size="20" tabindex="5"/><br/>
    <input type="submit" name="register" id="register"
	   value="Register" tabindex="6"/>
  </fieldset>
</form>
