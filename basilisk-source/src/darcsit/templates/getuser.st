<script type="text/javascript">
/* <![CDATA[ */
  \$.ajaxSetup({ cache: false });
  \$.get("$base$/_user", {}, function(username, status) {
     \$("#logged_in_user").text(username);
     if (username == "") {  // nobody logged in
        \$("#logoutlink").hide();
        \$("#loginlink").show();
     } else {
        \$("#logoutlink").attr("href", "$base$/_logout?_=" + Date.now());
        \$("#logoutlink").show();
        \$("#loginlink").hide();
     };
   });
/* ]]> */
</script>
