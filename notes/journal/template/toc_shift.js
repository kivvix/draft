#!/usr/bin/env node

var fs = require("fs");
const jsdom = require("jsdom");
const { JSDOM } = jsdom;

const argv = require('yargs')
  .scriptName("toc_shift")
  .usage("$0 -i inputfile -o outputfile")
  .option('input',{
    alias: 'i',
    describe: 'input file where need to move `#TOC` outside of `<main>` and add `.markdown-body` class to `<main>`'
  })
  .option('output',{
    alias: 'o',
    describe: 'output file with all modifications'
  })
  .demandOption(['input','output'],"Please provide both input and output arguments to work this tool")
  .help()
  .alias('help','h')
  .argv

JSDOM.fromFile( argv['input'] ,{}).then(dom => {
  // extract #TOC
  let toc  = dom.window.document.querySelector("#TOC");
  let main = dom.window.document.querySelector("main");

  // add checkbox and put #TOC in label
  let cb_TOC = dom.window.document.createElement('input');
  cb_TOC.type = 'checkbox';
  let cb_TOC_id = 'cb_TOC'
  cb_TOC.name = cb_TOC_id; cb_TOC.id= cb_TOC_id;

  let label = dom.window.document.createElement('label');
  label.htmlFor = cb_TOC_id;
  label.appendChild(toc)

  main.parentNode.insertBefore(label,main);
  main.parentNode.insertBefore(cb_TOC,label);

  // add markdown style to main
  main.classList.add("markdown-body");

  // write file
  fs.writeFile(argv['output'], dom.serialize(), (err) => {
    if (err) console.log(err);
    console.log("Successfully TOC shifting!");
  });

});
