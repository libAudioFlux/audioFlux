function setGithubUrl() {
    var oGitA = document.getElementsByClassName('fa fa-github')[0];
    oGitA.href = 'https://github.com/libAudioFlux/audioFlux';
    oGitA.innerHTML = ' Github';
    oGitA.target = '_blank';
}

window.onload = () => {
    setGithubUrl()
}