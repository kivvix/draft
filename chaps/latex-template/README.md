# Mes chapitres de thèse

> Il s'agit d'une reprise du template LaTeX de l'ED-MathSTIC dont le dépôt est sur [GitLab](https://gitlab.inria.fr/ed-mathstic/latex-template) (étrangement c'est sur celui de l'INRIA plutôt que celui de l'Université de Rennes 1), que j'ai un peu adapté à ma sauce.

## Structure du dépôt

* `main.tex` contient le squelette du document, aucun texte du manuscrit n'est présent dans ce fichier
* `these-dbl.cls` contient les dépendances, les paramètres de la bibliographie dont le style de citation et les paramètres de mise en page globale du manuscrit et plus particulièrement des deux couvertures
* `_couv/pagedegarde.tex` contient les variables à remplir par l'auteur pour compléter la page de garde, ces variables sont utilisées par `\maketitle` redéfini dans `these-dbl.cls`
* `couv/resume.tex` contient les variables à remplir par l'auteur pour compléter la quatrième de couverture, ces variables sont utilisées par les macros définies dans `these-dbl.cls`

Le Makefile vous aide à compiler le latex et la bibliographie en un pdf (détails plus bas). J'ai ajouté également un `make min` pour utiliser `ghostscript` pour minimiser le poids du PDF (pas très efficace sur ce fichier, la compression pourrait être plus agressive si besoin).
Les autres dossiers contiennent chacun un chapitre du document, ainsi :

* `chap1` contient(dra) l'article avec Lukas et Nicolas, avec peut-être (si j'ai le temps) quelques détails supplémentaires sur le code Python qui m'a permi d'obtenir automatiquement des estimations numériques de CFL.
* `chap2` contient la partie 1dx-1dv déjà écrite il y a environ 1 an
* `chap3` contient les débuts du 1dz-3dv, tout l'article n'a pas encore été repris.
* `intro` contiendra l'introduction
* `conclu` contiendra une possible conclusion

## Mes petits ajouts

En plus d'un simple `make min` j'ai aussi ajouté une variable `\localPath` dans le fichier principal de chaque chapitre. Ainsi l'inclusion d'une figure (quelque soit le chapitre) se fait via :

```
  \includegraphics{\localPath/figures/mafigure.png}
```

> Le contenu des dossiers `*/figures` est vérifié par le `Makefile` pour savoir si une nouvelle compilation est nécessaire, il est dont important que les figures soient présentes dans ces dossiers.

J'ai également ajouté la possibilté d'avoir des appendices propres à un chapitre (pour regrouper dans le chapitre des annexes qui lui sont associées). Cela est utile pour l'instant que dans le chapitre 2.

La bibliographie est commune à tous les chapitres, mais celle-ci indique les pages où sont référées les articles cités. J'ai mis tout mon fichier de biblio depuis le début de ma thèse, donc il y a un peu de ménage à y faire (prenons pour exemple des pdf de travail).

