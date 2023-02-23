# Hello genome paper collaborators!! 

This README.md file is sitting in a dedicated directory on Cedar, for committing files to GitHub related to our creek chub genome manuscript. Let me know if you have a better setup design in mind as I'm not all that familiar with Github but hopefully this will work well!

Directory path: 
> /project/rrg-emandevi/creekchub_genome/

## Commits

My idea is that we use a hard link for any files we want to add to the repo, so that we can work on them in our own directories or wherever, but have one dedicated directory from which we commit files. 
I used `cp -l /path/to/file .`to make a hard link, while in the `creekchub_genome` directory.

These lines of code should work for making commits and you can commit as many files as you'd like at a time:

```
git add <filename(s)> 
git commit -m <message>
git push 
```

## SSH Key

Be aware that you may have to make an ssh key in your GitHub profile before you're able to commit things. To do so, type:
`ssh-keygen -t rsa -C "your_github_email@email.com"`

Then enter through the options to add the key to the default location and not add a passcode. Then type:
`cat /home/username/.ssh/id\_rsa.pub`(if this was the location that it suggested to add the key to)

Then copy the key that it prints and paste it into the key box of the "SSH key" page, under profile settings.
