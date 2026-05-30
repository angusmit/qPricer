/ core/cfg.q - configuration loader (v0.57, ARCHITECTURE.md section 2)
/ ----------------------------------------------------------------------------
/ Populates the .cfg namespace. Loads config/base.q (defaults = the values
/ previously hardcoded in the modules), then optionally config/{env}.q to
/ override a subset, where the environment is selected by the QPRICER_ENV
/ variable. When QPRICER_ENV is unset (the default / test path), only base is
/ loaded, so behaviour is byte-identical to before this layer existed:
/   base.q + no override == prior hardcoded behaviour, exactly.
/ Loaded FIRST by core/init.q so .cfg is populated before any module that reads
/ it. Modules read .cfg lazily inside their default-config functions, so even a
/ load-time call would see a populated .cfg.
/ ----------------------------------------------------------------------------
.cfg.env:`$getenv `QPRICER_ENV;
system "l config/base.q";
if[not null .cfg.env;
    @[{[envName] system "l config/",envName,".q"};
      string .cfg.env;
      {[err] -1 "cfg: no usable override file for QPRICER_ENV; using base only (",err,")";}]
    ];
