# Adding password to a Jetstream instance

To add a password to your Jetstream instance, 'Open the Web Shell' from the instance page.

```
sudo passwd [username]
```
Enter a password when prompted. The letters will not display when you type, so do not be alarmed.

![](images/password-change.png)

Exit out of the Web Shell.

Open your terminal and login:

```
ssh [username]@[ip.address] 
```

Type your password when prompted.

![](images/jetstream_login.png)
