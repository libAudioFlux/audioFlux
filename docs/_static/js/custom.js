function setGithubUrl() {
    var oGitA = document.getElementsByClassName('fa fa-github')[0];
    oGitA.href = 'https://github.com/libAudioFlux/audioFlux';
    oGitA.innerHTML = ' Github';
    oGitA.target = '_blank';
}

function setVersion() {
    fetch(DOCUMENTATION_OPTIONS.switcher_version_json_url)
        .then((res) => {
            return res.json();
        })
        .then((data) => {
            const dl = document.getElementById("rst-version-dl");
            data.forEach((entry) => {
                if (!("name" in entry)) {
                    entry.name = entry.version;
                }
                const node = document.createElement("a");
                node.textContent = `${entry.name}`;
                node.setAttribute("href", `${entry.url}`)

                const dd = document.createElement('dd');
                dd.appendChild(node);
                dl.appendChild(dd);
            })
        })
}

window.onload = () => {
    setGithubUrl();
    setVersion();
}