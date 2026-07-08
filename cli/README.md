# Overview

CLI for running `bpact` with the option for posterior inference using:
- exact integration using a pseudo likelihood
- resampling (MCMC) using an exact likelihood

## TL;DR

```bash
pixi install
cp configs/template.yaml configs/analysis.yaml
# edit configs/analysis.yaml
pixi run bpact --config configs/analysis.yaml
```

## Environment

The `bpact` CLI uses [`pixi`](https://pixi.prefix.dev/latest/#installation) to handle the environment. [Download `pixi`](https://pixi.prefix.dev/latest/#installation) then run from `bpact/cli`:

```bash
pixi install
```

## Running

> [!NOTE]
> The `bpact` CLI requries GWAS summary statistics


