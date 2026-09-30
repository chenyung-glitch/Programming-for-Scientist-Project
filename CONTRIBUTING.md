# 🚀 Team Development Guide

Welcome! Please follow these guidelines to keep our codebase clean and prevent broken builds.

---

## 📥 1. Project Setup (Cloning to Local Directory)

Before you can work, you need to copy the remote repository onto your computer. Choose **one** of the methods below:

### Option A: Using the Terminal (Recommended)

1. Open your terminal or command prompt.
2. Navigate to the folder where you want to keep your school projects (e.g., `cd Documents/CS-Projects`).
3. Run the clone command using our repository URL:
   ```bash
   git clone https://github.com/chenyung-glitch/Programming-for-Scientist-Project.git
   or
   git clone git@github.com:chenyung-glitch/Programming-for-Scientist-Project.git # if you have setup your ssh key
   ```
4. Move into the newly created project folder:
   ```bash
   cd Programming-for-Scientist_Project
   ```

### Option B: Using the GitHub Desktop App

1. Open **GitHub Desktop** and sign into your account.
2. Click on **File** (top menu) -> **Clone repository...**
3. Select our repository from the list (or click the **URL** tab and paste our repository link).
4. Under **Local Path**, click *Choose...* to select the exact folder on your computer where you want the project saved.
5. Click **Clone**.

---

## 🌿 2. Git Flow Procedure

### 1. Starting a Task

Always create a new branch for your work on your local copy. Never work directly on `main`.

```bash
git checkout main # go to the main branch on your local copy
git pull origin main # pull the most recent version of the project from the remote github
git checkout -b feature/your-feature-name # create and go to a new branch based on what you are working on
```

### 2. Modifying & Pushing Files

Keep your commits small, focused, and give them clear messages.

```bash
git status # check all the files that have been modified and added before you add to your queue
git checkout feature/your-feature-name # move to the branch where you want the changes to be at
git add . # add all changes that you have made to your queue
git commit -m "feat: description of changes" # save all your changes in the queue with a message 
# you can keep using git add and git commit to save your progress everytime you do a small task
git push origin feature/your-feature-name # push your changes to the remote github to be approved
hello
```

### 3. Merging

1. Go to GitHub and open a **Pull Request (PR)**.
2. Tag at least **one teammate** for a code review.
3. Once approved, merge the PR and clean up your local machine:

```bash
git log --oneline # checks the branches and change history of the repo
git checkout main # move to main
git pull origin main # pull the merged version
git branch -d feature/your-feature-name # delete the previous branch on your local machine because it has now been merged
```

---

## 🐍 3. Managing Python Packages

Our virtual environment folder (`.venv/`) is ignored by Git. Instead, we track packages using `requirements.txt`.
I have setup the requirements.txt, so all you need to do to setup your own local python environment as below.

### Setting Up Python Environment

```bash
git pull origin main # pull the latest code
python -m venv .venv # create their own local environment
# or python3 -m venv .venv if above doesn't work
source .venv/bin/activate # activate it
# or .venv\Scripts\activate on Windows
pip install -r requirements.txt # install everything from the shared recipe list
```

### Installing a New Package

If you need a new dependency (e.g., `requests`), follow these steps on your feature branch:

1. Ensure your environment is active: `source .venv/bin/activate` (or `.venv\Scripts\activate` on Windows).
2. Install the package: `pip install requests`
3. Update the recipe file: `pip freeze > requirements.txt`
4. Commit both your code changes and the updated `requirements.txt`.

### Updating Your Local Environment

When someone else adds a package and you pull their changes down, sync your environment by running:

```bash
pip install -r requirements.txt
```

---

## 4. Adding to README.md

Make sure that after working on a file and doing all the add and commits needed, add to the README.md with a brief description of the file before pushing to main!
